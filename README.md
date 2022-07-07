### NPMCV

A Python script for segmenting cells in microscope images and calculating the coefficient of variation (CV).

## Installation

Windows
```
$ python -m build
$ pip install
```
Mac
```
$ python setup.py build
$ python3 setup.py install 
```
## How to use

```
Usage:
    npmcv [<directory>...]
```
where <directory> is a one or more folders containing `.lif` images.

```
Outputs:
   <directory>.csv:     output file
   <directory>_RAW.csv: output data without outliners removed  
   <directory>_OUT.csv: the outliners removed from data
```
Columns in `RAW` output corresponds to the CV values calculated for each segmented cell in an individual image.
Columns in standard output file are grouped by `.lif` file 

# Additional Outputs:
`<directory>/dapi_seg/<img>_seg.png:`

Generates an image overlaying the cell segmentation on the original image
to validate the segmentation.

`<directory>/inv_cells/<img>_cell<n>.png:`

The actual cells where CV was calculated, `<n>` corresponds to the row number in the csv files.
