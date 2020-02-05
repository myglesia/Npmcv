

## Installation

change cwd to `Npmcv`
```
$ python setup.py build
$ python setup.py install
```
## How to use

```
Usage:
    npmcv [<directory>...]
```
where <directory> is a one or more folders containing tif images. 

```
Output:
   <directory>.csv:     output file
   <directory>_RAW.csv: output data without outliners removed  
   <directory>_OUT.csv: the outliners removed from data
```
