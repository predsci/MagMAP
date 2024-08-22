<img width=500 src="doc/magmap_logo.png" alt="HipFT" />  

# MagMAP: Magnetic Mapping And Processing
    
[Predictive Science Inc.](https://www.predsci.com)  
 
--------------------------------  

## Overview  
  
MagMAP is a python package for downloading, mapping, and processing solar magnetograms.  
It is an integral component 
of the Open-source Flux Transport [OFT](https://github.com/predsci/OFT) model.  

## Requirements

This package requires specific versions of python packages.
To facilitate this, we have provided a conda recipe file to 
create a python environment that MagMAP can run in.
Some packages/versions are in the conda-forge channel. 
To add this channel, first check if it is there:  
`conda config --show channels`  
If it is not there, add conda-forge:  
`conda config --append channels conda-forge`  
The MagMAP conda environment can then be built using:  
`conda env create -f rsrc/conda_recipe_magmap.yml`  
To activate the environment, use:  
`conda activate mapmag`  
  

## Getting Started
### Data Download
See the following files for examples of how to download a set of HMI line-of-sight
magnetorgrams - magmap/examples/Download_HMI.py, and then update the directory with
the most recent magnetograms - magmap/data/scripts/Update_local_HMI.py.  The primary 
function of this code is to specify a time range and cadence, download 
available data that matches the specification, and keep everything in an orderly
directory structure with updated index file.  The update script (Update_local_HMI.py)
writes an index file all-files.csv that lists the timestamp and relative path of 
each data file.

### Mapping Data
Once a local directory of HMI data exists, the script magmap/examples/Map_HMI_dir.py
can be used to generate maps.  Additionally, the script 
magmap/maps/scripts/Update_Map_HMI.py can be used to update maps as new data comes 
in.  The mapping process has five primary steps
- Load the disk FITS file into a LosMagneto object
- High resolution interpolation from disk to map coordinates
- Reduce resolution by integration
- Set assimilation weights for HipFT
- Save in an HDF5 file format

## Contact

James Turtle ([jturtle@predsci.com](mailto:jturtle@predsci.com))

## Contributors

- James Turtle
- Ronald M. Caplan
- Jon A. Linker
- Cooper Downs
