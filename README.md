<img height=40 src='http://my5C.umassmed.edu/images/3DG.png' title='3D-Genome' />
&nbsp;&nbsp;
<img height=30 src='http://my5C.umassmed.edu/images/dekkerlabbioinformatics.gif' />
&nbsp;&nbsp;
<img height=30 src='http://my5C.umassmed.edu/images/umasslogo.gif' />

# c-world-encode

h5 helper scripts for encode/dekkerlab hi-c data.

convert hdf5 (h5) files into tsv (tab seperated value) text (matrix) files.

```
c-world-encode/
	hdf2tab.py - Extract c-data from HDF5 file into TXT (matrix.gz)
```

## Full Documentation

See the Wiki for full documentation of encode/dekkerlab Hi-C HDF5 file format.
[c-world-encode wiki](https://github.com/blajoie/c-world-encode/wiki)

[HDF5 spec](https://github.com/blajoie/c-world-encode/wiki/H5-Spec)

[HDF5 helper script (python), hdf2tab.py](https://github.com/blajoie/c-world-encode)
numpy/scipy/h5py required. Python 2.7+

[hdf2tab.py usage](https://github.com/blajoie/c-world-encode#usage</a>)

Hi-C data is pushed through the dekkerlab standard pipeline (to be available on GIT soon) and binned at multiple bin sizes (resolutions)
10mb
2.5mb
1mb
500kb
250kb
100kb
40kb

See the [Wiki](https://github.com/blajoie/c-world-encode/wiki) for full documentation, examples, operational details and other information.

## Communication

- [Bryan Lajoie](https://github.com/blajoie)
- [Noam Kaplan](https://github.com/NoamKaplan)
- Twitter: [@my5C](https://twitter.com/my5C)

## What does it do?

hdf2tab can read/subset a hdf5 file containing Hi-C data from the encode/dekkerlab and convert it into a my5C fomatted tsv matrix file.

## Usage

```
$ git clone git@github.com:blajoie/c-world-encode.git
$ cd c-world-encode/
$ python hdf2tab.py

usage: hdf2tab.py [-h] -i INFILE [-info] [-v] [-o OUTFILE] [-cis]
                  [-chrs SELECTED_CHRS [SELECTED_CHRS ...]]
                  [-z ZOOM_COORDS [ZOOM_COORDS ...]] [-m BLOCKMEM]
                  [-p PRECISION] [-cl CHRLISTFILE] [-bl BINLISTFILE]

Extract c-data from HDF5 file into TXT (matrix.gz)

optional arguments:
  -h, --help            show this help message and exit
  -i INFILE, --input INFILE
                        interaction matrix hdf5 file (default: None)
  -info                 interaction matrix hdf5 file (default: False)
  -v, --verbose         Increase verbosity (specify multiple times for more)
                        (default: None)
  -o OUTFILE, --output OUTFILE
                        interaction matrix output file (default: None)
  -cis                  extract cis maps only (default: False)
  -chrs SELECTED_CHRS [SELECTED_CHRS ...]
                        subset of chromosomes to extract (default: *)
  -z ZOOM_COORDS [ZOOM_COORDS ...], --zoom ZOOM_COORDS [ZOOM_COORDS ...]
                        zoom coordinate (can only select symmetrical subsets)
                        (default: None)
  -m BLOCKMEM, --bmem BLOCKMEM
                        block size for extracting (default=hdf chunk size)
                        (default: None)
  -p PRECISION          output precision (# of digits) (default: 4)
  -cl CHRLISTFILE, --chrlist CHRLISTFILE
                        chromosome list output file (default: None)
  -bl BINLISTFILE, --binlist BINLISTFILE
                        bin position output file (default: None)
usage: hdf2tab.py [-h] -i INFILE [-info] [-v] [-o OUTFILE] [-cis]
                  [-chrs SELECTED_CHRS [SELECTED_CHRS ...]]
                  [-z ZOOM_COORDS [ZOOM_COORDS ...]] [-m BLOCKMEM]
                  [-p PRECISION] [-cl CHRLISTFILE] [-bl BINLISTFILE]
hdf2tab.py: error: argument -i/--input is required

```

## Bugs and Feedback

For bugs, questions and discussions please use the [Github Issues](https://github.com/blajoie/c-world-encode/issues).

## LICENSE

Licensed under the Apache License, Version 2.0 (the 'License');
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

<http://www.apache.org/licenses/LICENSE-2.0>

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an 'AS IS' BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

