
from __future__ import print_function

import numpy as np
import scipy as sp
import pdb
import h5py
import sys
import argparse
import logging
import time
import gzip
import re
import os
import math
from psutil import virtual_memory
from collections import defaultdict

def main():

    parser=argparse.ArgumentParser(description='Extract c-data from HDF5 file into TXT (matrix.gz)',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--input', dest='infile', type=str, required=True, help='interaction matrix hdf5 file')
    parser.add_argument('-info', dest='info', action='store_true', help='interaction matrix hdf5 file')
    parser.add_argument('-or','--output_relative', dest='output_relative', action='store_true', help='output file relative to input file path')
    parser.add_argument('-v', '--verbose', dest='verbose',  action='count', help='Increase verbosity (specify multiple times for more)')
    parser.add_argument('-o', '--output', dest='outfile', type=str, help='interaction matrix output file')
    parser.add_argument('-cis', dest='cis_mode', action='store_true', help='extract cis maps only')
    parser.add_argument('-chrs', dest='selected_chrs', nargs='+', type=str, default='*', help='subset of chromosomes to extract')
    parser.add_argument('-maxdim', dest='max_dimension', type=int, default=300000, help='maximum dimension of allxall matrix - else cis only')
    parser.add_argument('-z', '--zoom', dest='zoom_coords', nargs='+', type=str, help='zoom coordinate (can only select symmetrical subsets)')
    parser.add_argument('-m','--bmem', dest='blockmem', type=int, help='block size for extracting (default=hdf chunk size)')
    parser.add_argument('-p', dest='precision', type=int, default=4, help='output precision (# of digits)')
    parser.add_argument('-cl', '--chrlist', dest='chrlistfile', type=str, help='chromosome list output file')
    parser.add_argument('-bl', '--binlist', dest='binlistfile', type=str, help='bin position output file')
    
    #parser.print_help()
    #usage = "usage: %prog [options] arg1 arg2"
    
    args=parser.parse_args()

    infile=args.infile
    info=args.info
    output_relative=args.output_relative
    verbose=args.verbose
    outfile=args.outfile
    cis_mode=args.cis_mode
    selected_chrs=args.selected_chrs
    max_dimension=args.max_dimension
    zoom_coords=args.zoom_coords
    blockmem=args.blockmem
    precision=args.precision
    chrlistfile=args.chrlistfile
    binlistfile=args.binlistfile
    
    log_level = logging.WARNING
    if verbose == 1:
        log_level = logging.INFO
    elif verbose >= 2:
        log_level = logging.DEBUG
    logging.basicConfig(level=log_level)
    
    verbose = info if info else verbose
    verboseprint = print if verbose else lambda *a, **k: None
    format_func=("{:."+str(precision)+"f}").format
    
    verboseprint("\n",end="")
    
    infile_name=os.path.basename(infile)
    if output_relative:
       infile_name=infile
    inhdf=h5py.File(infile,'r')
    
    # attrs
    genome=inhdf.attrs['genome'][:]
    # datasets
    bin_positions=inhdf['bin_positions'][:]
    chr_bin_range=inhdf['chr_bin_range'][:]
    chrs=inhdf['chrs'][:]
    # matrix shape
    nrow=inhdf['interactions'].shape[0]
    ncol=inhdf['interactions'].shape[1]
    
    # calculate optimal block size
    itx_dtype=inhdf['interactions'].dtype
    itx_dtype_size=itx_dtype.itemsize
    hdf_blocksize=inhdf['interactions'].chunks[0]
    hdf_blockmem=byte_to_megabyte(hdf_blocksize*(ncol*itx_dtype_size))
    
    blocksize,blockmem=get_block_size(itx_dtype_size,hdf_blocksize,hdf_blockmem,blockmem,[nrow,ncol])    
    verboseprint("blocksize",blocksize)
    verboseprint("blockmem",blockmem)
    verboseprint("precision",precision)
    verboseprint("")
    
    # ensure symmetrical
    if nrow!=ncol:
        sys.exit('error: non-symmetrical matrix found!')
    n=nrow=ncol
    
    # build chr lookup dict
    chr_dict={}
    for i,c in enumerate(chrs):
        chr_dict[c]=i
    
    # dump hdf info if selected
    if(info):
        verboseprint("inputFile",infile,sep="\t")
        verboseprint("inputFileName",infile_name,sep="\t")
        verboseprint("matrix shape\t",nrow," x ",ncol,sep="")
        verboseprint("assembly",genome,sep="\t")
        verboseprint("h5 chunk",inhdf['interactions'].chunks[0],sep="\t")
        verboseprint("user chunk",blocksize,sep="\t")
        
        verboseprint("\nchrs",sep="\t")
        for i,c in enumerate(chrs):
            cbr=chr_bin_range[chr_dict[c]]
            start,end=bin_positions[cbr[0]][1],bin_positions[cbr[1]][2]
            size=(end-start)+1
            nbins=(cbr[1]-cbr[0])+1
            verboseprint("\t",i,"\t",c,":",start,"-",end,"\t(",size,")\t",cbr,"\t",nbins,sep="")
            
        verboseprint("")
        quit()
    
    # warn user that (txt) matrix > 300,000 row/col is _excessive_
    if(n>max_dimension):
        logging.warning("large matrix! %d > %d.\n\tenforcing cis only mode!\n" % (n,max_dimension))
        cis_mode=1
    
    if(outfile==None):
        outfile=re.sub(".hdf5", "", infile_name)
    outfile=re.sub(".matrix", "", outfile)
    outfile=re.sub(".gz", "", outfile)

    # process zoom coordinates
    verboseprint("zoom coordinates")
    zoom_chrs=list()
    zoom_dict=defaultdict(list)
    if(zoom_coords!=None):
        for z in zoom_coords:
            zoom_coord=split_zoom_coord(z)
            
            if zoom_coord==None:
                verboseprint("\t",z," *invalid*",sep="")
                continue
                
            zoom_chr,zoom_start,zoom_end=zoom_coord
            verboseprint("\t",zoom_chr,":",zoom_start,"-",zoom_end,sep="")
            zoom_chrs += [zoom_chr]
            zoom_dict[zoom_chr].append(zoom_coord)
        verboseprint("\n")
    
    # process selected chromosomes
    verboseprint("selected chromosomes")
    if selected_chrs!="*":        
        selected_chrs=de_dupe_list(selected_chrs+zoom_chrs)
        for c in selected_chrs:
            if not c in chr_dict:
                sys.exit('specificed chr '+c+'not found in file!')
        selected_chrs=sorted(selected_chrs,key=lambda x:chr_dict[x]) # ensure selected_chrs are sorted same as the HDF
        for c in selected_chrs:
            verboseprint("\t",c,"\n",end="")
    else:
        selected_chrs=chrs
        verboseprint("\t*\n",end="")
        
    verboseprint("")
    
    bin_mask=np.zeros(n,dtype=bool)
    
    # dump hdf, chr x chr (multiple matrix)
    
    if(cis_mode == 1):
    
        verboseprint("cis only mode\n")
        
        for c in selected_chrs:
            c_ind=chr_dict[c]
            r=chr_bin_range[chr_dict[c]]
            
            # reset bin_mask to all zeros
            bin_mask=np.zeros(n,dtype=bool)
            if c in zoom_dict:
                zoom_coord_arr=zoom_dict[c]
                for zoom_coord in zoom_coord_arr:
                    tmp_bin_positions=bin_positions[r[0]:r[1]+1]
                    for i,b in enumerate(tmp_bin_positions):
                        if b[2] < zoom_coord[1]: continue
                        if b[1] > zoom_coord[2]: break
                        overlap=is_overlap([zoom_coord[1],zoom_coord[2]], [b[1],b[2]])
                        if(overlap > 0):
                            bin_mask[r[0]+i]=True
                    n2=np.sum(bin_mask[r[0]:r[1]+1])
                    verboseprint("\t",c," zoom subset ",n2,"x",n2,sep="")
            else:
                bin_mask[r[0]:r[1]+1]=True
                n2=np.sum(bin_mask[r[0]:r[1]+1])
                verboseprint("\t",c," chr subset ",n2,"x",n2,sep="")
                
            # interaction matrix output
            n2=np.sum(bin_mask)
            
            if n2 > max_dimension:
                logging.warning("sub-matrix too large! [%s] %d > %d.\n" % (c,n,max_dimension))
                continue
                
            if not verbose: print(c," - writing (",n2,"x",n2,") matrix",sep="")
            m_out_fh=gzip.open(outfile+'__'+c+'.matrix.gz',"wb")

            header=[str(i)+'|'+genome+'|'+str(chrs[bin_positions[i,0]])+':'+str(bin_positions[i,1])+'-'+str(bin_positions[i,2]) for i in np.nonzero(bin_mask)[0]]
            print(str(n2)+"x"+str(n2)+"\t"+"\t".join(header),file=m_out_fh)
                       
            k=0
            tmp_chr_bin_mask=bin_mask[r[0]:r[1]+1]
            c_end=max(np.nonzero(tmp_chr_bin_mask))[-1]
            c_start=max(np.nonzero(tmp_chr_bin_mask))[0]
            tmp_chr_bin_mask=tmp_chr_bin_mask[c_start:c_end+1]
            c_start += r[0]
            c_end += r[0]
                        
            for i in xrange(c_start,c_end+1,blocksize):
                b=min(c_end+1-i,blocksize)
                tmp_bin_mask=tmp_chr_bin_mask[i-c_start:i-c_start+b]
                    
                verboseprint("\r",""*20,"\r\tloading block (",i,":",i+b,") ... ",sep="",end="\r")
                if verbose: sys.stdout.flush()
                
                current_block=inhdf['interactions'][i:i+b,:][:,bin_mask][tmp_bin_mask,:]
                
                for j in xrange(current_block.shape[0]):
                    print(header[k]+"\t"+"\t".join(map(format_func,current_block[j,:])),file=m_out_fh)
                    m_out_fh.flush()
                    
                    pc=((float(k)/float((n2)))*100)
                    verboseprint("\r\t"+str(k)+" / "+str(n2)+" ["+str("{0:.2f}".format(pc))+"%] complete ... ",end="\r")
                    if verbose: sys.stdout.flush()
                    
                    k+=1

            m_out_fh.close()
            
            verboseprint('\r',end="")
            pc=((float(n2)/float((n2)))*100)
            verboseprint("\t"+str(n2)+" / "+str(n2)+" ["+str("{0:.2f}".format(pc))+"%] complete",end="")
            if verbose: sys.stdout.flush()
            
            verboseprint("")
            
            if (chrlistfile!=None):
                verboseprint("writing chr list file")
                c_out_fh=open(chrlistfile+'__'+c+'.cl',"w")
                for c in selected_chrs:
                    print(c,file=c_out_fh)
                c_out_fh.close()
                verboseprint("\tdone\n")

            if (binlistfile!=None):
                verboseprint("writing bin list file")
                b_out_fh=open(binlistfile+'__'+c+'.bp',"w")
                for i in np.nonzero(bin_mask)[0]:
                    print(str(i)+"\t"+chrs[bin_positions[i,0]]+"\t"+str(bin_positions[i,1])+"\t"+str(bin_positions[i,2]),file=b_out_fh)
                b_out_fh.close()
                verboseprint("\tdone\n")
                
            verboseprint("")
    else:
        
        # dump hdf, all x all (one matrix)
        
        verboseprint("all mode\n")
        
        for c in selected_chrs:
            r=chr_bin_range[chr_dict[c]]
            if c in zoom_dict:
                zoom_coord_arr=zoom_dict[c]
                for zoom_coord in zoom_coord_arr:
                    tmp_bin_positions=bin_positions[r[0]:r[1]+1]
                    for i,b in enumerate(tmp_bin_positions):
                        if b[2] < zoom_coord[1]: continue
                        if b[1] > zoom_coord[2]: break
                        overlap=is_overlap([zoom_coord[1],zoom_coord[2]], [b[1],b[2]])
                        if(overlap > 0):
                            bin_mask[r[0]+i]=True
                    n2=np.sum(bin_mask[r[0]:r[1]+1])
                    verboseprint("\t",c," zoom subset ",n2,"x",n2,sep="")
            else:
                bin_mask[r[0]:r[1]+1]=True
                n2=np.sum(bin_mask[r[0]:r[1]+1])
                verboseprint("\t",c," chr subset ",n2,"x",n2,sep="")
        
        verboseprint("")
        
        # interaction matrix output
        n2=np.sum(bin_mask)
        
        if n2 > max_dimension:
            
            logging.warning("matrix too large! %d > %d.\n" % (n,max_dimension))
            
        else:
            
            if not verbose: print("all - writing (",n2,"x",n2,") matrix",sep="")
            m_out_fh=gzip.open(outfile+'.matrix.gz',"wb")
            
            # create header list
            header=list()
            for c in selected_chrs:
                tmp_bin_mask=np.zeros(n,dtype=bool)
                r=chr_bin_range[chr_dict[c]]
                tmp_chr_bin_mask=bin_mask[r[0]:r[1]+1]
                tmp_bin_mask[r[0]:r[1]+1]=tmp_chr_bin_mask
                tmp_header=[str(i)+'|'+genome+'|'+str(chrs[bin_positions[i,0]])+':'+str(bin_positions[i,1])+'-'+str(bin_positions[i,2]) for i in np.nonzero(tmp_bin_mask)[0]]
                header += tmp_header
                
            print(str(n2)+"x"+str(n2)+"\t"+"\t".join(header),file=m_out_fh)
          
            k=0
            for c in selected_chrs:
                r=chr_bin_range[chr_dict[c]]
                tmp_chr_bin_mask=bin_mask[r[0]:r[1]+1]
                c_end=max(np.nonzero(tmp_chr_bin_mask))[-1]
                c_start=max(np.nonzero(tmp_chr_bin_mask))[0]
                tmp_chr_bin_mask=tmp_chr_bin_mask[c_start:c_end+1]
                c_start += r[0]
                c_end += r[0]

                for i in xrange(c_start,c_end+1,blocksize):
                    b=min(c_end+1-i,blocksize)
                    tmp_bin_mask=tmp_chr_bin_mask[i-c_start:i-c_start+b]
                        
                    verboseprint("\r",""*20,"\r\tloading block (",i,":",i+b,") ... ",sep="",end="\r")
                    if verbose: sys.stdout.flush()
                    
                    current_block=inhdf['interactions'][i:i+b,:][:,bin_mask][tmp_bin_mask,:]
                    
                    for j in xrange(current_block.shape[0]):
                        print(header[k]+"\t"+"\t".join(map(format_func,current_block[j,:])),file=m_out_fh)
                        m_out_fh.flush()
                        
                        pc=((float(k)/float((n2)))*100)
                        verboseprint("\r\t"+str(k)+" / "+str(n2)+" ["+str("{0:.2f}".format(pc))+"%] complete ... ",end="\r")
                        if verbose: sys.stdout.flush()
                        
                        k+=1
                
            m_out_fh.close()
            
            verboseprint('\r',end="")
            pc=((float(n2)/float((n2)))*100)
            verboseprint("\t"+str(n2)+" / "+str(n2)+" ["+str("{0:.2f}".format(pc))+"%] complete",end="")
            if verbose: sys.stdout.flush()
            
            verboseprint("")
            
            if (chrlistfile!=None):
                verboseprint("writing chr list file")
                c_out_fh=open(chrlistfile+'.cl',"w")
                for c in selected_chrs:
                    print(c,file=c_out_fh)
                c_out_fh.close()
                verboseprint("\tdone\n")

            if (binlistfile!=None):
                verboseprint("writing bin list file")
                b_out_fh=open(binlistfile+'.bp',"w")
                for i in np.nonzero(bin_mask)[0]:
                    print(str(i)+"\t"+chrs[bin_positions[i,0]]+"\t"+str(bin_positions[i,1])+"\t"+str(bin_positions[i,2]),file=b_out_fh)
                b_out_fh.close()
                verboseprint("\tdone\n")
                
            verboseprint("")

def split_zoom_coord(z):
    """validate and split zoom coordinate.
    coordinates must be UCSC formatted.
    e.g. chr1:500-1000
    chr(colon)start(hyphen)end where start <= end
    """
    zoom_coord=re.search(r'(\S+):(\d+)-(\d+)',z)
    
    if zoom_coord==None:
        return None
        
    zoom_chr,zoom_start,zoom_end=zoom_coord.groups()
    zoom_start=int(zoom_start)
    zoom_end=int(zoom_end)
    
    if(zoom_start > zoom_end):
        return None
        
    return [zoom_chr,zoom_start,zoom_end]
        
def get_block_size(itx_dtype_size,hdf_blocksize,hdf_blockmem,blockmem,hdf_shape):
    """choose optimal hdf block size, based on user blockmem request.
    pixel dtype * ncol = mem of single row.
    Adjust block size (num of rows) to acheive blockmem mem usage
    """
    
    mem = virtual_memory()
    system_mem=byte_to_megabyte(mem.total)
    
    nrow,ncol=hdf_shape
    
    if blockmem==None:
        blocksize=hdf_blocksize
        blockmem=byte_to_megabyte(blocksize*(ncol*itx_dtype_size))
    else:
        blockmem = system_mem/2 if blockmem > system_mem else blockmem
        blockmem -= 256 # reserve 256MB for non-matrix items
        num_block_chunks=math.floor(blockmem/hdf_blockmem)
        if(num_block_chunks > 1):
            blocksize = int(hdf_blocksize*num_block_chunks)
        else:
            blocksize = hdf_blocksize
        blockmem=byte_to_megabyte(blocksize*(ncol*itx_dtype_size))
    return blocksize,blockmem

def de_dupe_list(input):
    """de-dupe a list, preserving order.
    """
    
    output = []
    for x in input:
        if x not in output:
            output.append(x)
    return output

def byte_to_megabyte(byte):
    """convert bytes into megabytes.
    """
    
    return round(((byte / 1000) / 1000),4) # megabyte
    # return round(((byte / 1024) / 1024),4) # mebibyte

def is_overlap(a, b):
    """test to for overlap between two intervals.
    """
    
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))

if __name__=="__main__":
    main()

   