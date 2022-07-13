### NPMCV

A Python script for segmenting cells in microscope images and calculating the coefficient of variation (CV).

## Installation


Download the latest `.whl` package from the release page and install with `pip`

```
$ pip install npmcv-2.1.3.whl
```

Alternatively you can build the package yourself from the source with `python -m build`



## How to Use

```
usage: npmcv [-h] [-D] [-n [NAME]] [-V] <path>

positional arguments:
  <path>                directory containing lif images

optional arguments:
  -h, --help            show this help message and exit
  -D, --no_imgs         Don't save overlays and individual cell images
  -n [NAME], --name [NAME]
                        Set basename for the output files
  -V, --version         show program's version number and exit
```


The script generates several output files
```
   NAME.csv:     cleaned output file
   NAME_RAW.csv: output data without outliners removed  
   NAME_OUT.csv: the outliners removed from data
```
By default the output files will be named after the input folder.
Columns in `RAW` output corresponds to the CV values calculated for each segmented cell in an individual image.
Columns in standard output file are grouped by `.lif` file 

# Additional Image Outputs:

By default the program will save processed cell images to verify segmentation.

`<directory>/dapi_seg/<img>_seg.png:`

Generates an image overlaying the cell segmentation on the original image
to validate the segmentation.

`<directory>/inv_cells/<img>_cell<n>.png:`

The actual cells where CV was calculated, `<n>` corresponds to the row number in the csv files.
