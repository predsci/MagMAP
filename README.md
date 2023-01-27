# OFTpy: Open source Flux Transport model

OFTpy is a python package for downloading and processing solar magnetograms.  In conjunction with
the Fortran codebase [HipFT](https://github.com/predsci/HipFT), it forms the Open-source Flux Transport
(OFT) package.  OFT is an open-source, GPU-enabled package that handles everything from data acquisition 
to simulation and analysis for solar flux transport models.

## Requirements

This package requires a custom python environment that includes Sunpy.  
At time of writing the version of Sunpy (3.1.6) available from the 
conda-defaults channel is fatally flawed.  So be sure that the conda-forge 
channel is added to your conda preferences before using the enclosed 
conda recipe file.  Check for conda-forge:
```
conda config --show channels
```
If needed, add conda-forge:
```
conda config --append channels conda-forge
```
Now we are ready to build the custom 'oft' conda environment:

```
conda env create --file conda_recipe_oft.yml
```

## Getting Started
### Data Download
See the following files for examples of how to download a set of HMI line-of-sight
magnetorgrams - oftpy/examples/Download_HMI.py, and then update the directory with
the most recent magnetograms - oftpy/data/scripts/Update_local_HMI.py.  The primary 
function of this code is to specify a time range and cadence, download 
available data that matches the specification, and keep everything in an orderly
directory structure with updated index file.  The update script (Update_local_HMI.py)
writes an index file all-files.csv that lists the timestamp and relative path of 
each data file.

### Mapping Data
Once a local directory of HMI data exists, the script oftpy/examples/Map_HMI_dir.py
can be used to generate maps for HipFT.  Additionally, the script 
oftpy/maps/scripts/Update_Map_HMI.py can be used to update maps as new data comes 
in.  The mapping process has five primary steps
- Load the disk FITS file into a LosMagneto object
- High resolution interpolation from disk to map coordinates
- Reduce resolution by integration
- Set assimilation weights for HipFT
- Save to HipFT specification in an HDF5 file format

## Contact

James Turtle ([jturtle@predsci.com](mailto:jturtle@predsci.com))

## Contributors

- James Turtle
- Ronald M. Caplan
- Jon A. Linker
- Cooper Downs
