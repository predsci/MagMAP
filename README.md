<img width=500 src="doc/magmap_logo.png" alt="HipFT" />  

# MagMAP: Magnetic Mapping And Processing
    
[Predictive Science Inc.](https://www.predsci.com)  
 
--------------------------------  

## OVERVIEW  
  
MagMAP is a python package for downloading, mapping, and processing solar 
magnetograms. It is an integral component 
of the Open-source Flux Transport [OFT](https://github.com/predsci/OFT) model.  

--------------------------------  
   
## HOW TO BUILD MAGMAP

First, download or clone this repository. For example:  
```
git clone https://github.com/predsci/MagMAP.git
```

MagMAP had been tested to work with specific versions of some python packages.  
In order to facilitate this requirement, a conda recipe file is provided in `rsrc/conda_recipe_magmap.yml`.  
See [here](https://docs.anaconda.com/miniconda) for instructions on installing a conda environment.  
Be sure to confirm conda is active in your current shell.  
As some packages/versions are in the conda-forge channel, first check that it is in your python environment with:  
```
conda config --show channels
```  
If `conda-forge` is not there, add it with:  
```
conda config --append channels conda-forge
```  
The MagMAP conda environment can then be built using:  
```
conda env create -f rsrc/conda_recipe_magmap.yml
```  
To activate the environment, use:  
```
conda activate magmap
```  
This environment should be active before each use of `magmap`.
  
While in the directory of the repository clone, install the `magmap` package with:
```  
pip install ${PWD}  
```

--------------------------------  

## HOW TO RUN MAGMAP

The standard use of MagMAP involves two steps:  data acquisition of magnetogram disk images, and mapping them into Carrington coordinate maps (along with converting the disk image data to radial magnetic field).  If using the conda magmap environment described above, ensure it is active before running `magmap`.

### Data Acquisition

The file `bin/magmap_get_data.py` is an executable script that downloads HMI line-of-site magnetograms. By default, given start and end timestamps will download HMI_M 720s images between those dates/times at a one hour cadence:
```
python magmap_get_data.py 2024-01-01T00:00:00 2024-01-02T00:00:00
```
The script will download the data and store them in an orderly directory structure with an index CSV file that lists the timestamp and relative path of each data file.  

The script has optional arguments to set the data cadence, search window, output directory, and index file name.  
Run `python magmap_get_data.py -h` for more information on these options.  
  
The script can also update an existing folder of disk data with the most recently available data by running the script with the output directory set to the existing folder.   
   
For another example script on how `magmap` can update an existing directory with the most recent HMI magnetograms, see `magmap/data/scripts/Update_local_HMI.py`.
  
### Processing and Mapping Data

Once a directory of HMI disk magnetogram image data exists, the script `bin/magmap_disk2map.py` is used to convert the LOS magnetogram images into radial magnetic field values (Br) and mapped to longitude-colatitude Carrington coordinate maps.  

The process has five primary steps:
  
 - Load the disk data and convert it into Br
 - Interpolate the Br disk data to a very high resolution "interpolation" Carrington coordinate map
 - Bin the high resolution map down to the chosen final map resolution (default 1024x512) using flux-preserving integration
 - Set a layer of default data assimilation weights for use with [HipFT](https://github.com/predsci/hipft) as well as an additional layer with mu=cos(t) where `t` is the disk-to-limb angle for use with custom assimilation functions in HipFT
 - Save the final three-layer maps in 3D HDF5 files

An example of running `magmap_disk2map.py`on a folder containing the output of `magmap_get_data.py` called here `magmap_data_disks` is:
```
python magmap_disk2map.py magmap_data_disks
```
The script will process and map the disk data and store them in an orderly directory structure with an index CSV file that lists the timestamp and relative path of each map file.  

The script has optional arguments to set the start and end date/time (for processing subsets), the theta and phi resolutions of the final maps, the assumed radius of the disk data (in solar radii), the theta and phi resolutions used for the high resolution interpolation map, output directory, and index file name.  Run `python magmap_disk2map.py -h` for more information on these options.

Once a map output folder has generated, it can be directly used with HipFT (see instructions there).

The script can also be used to update a pre-existing map output folder by rerunning it on the same disk data folder with the same output folder name.

For the purpose of automated data processing, an alternative script  `magmap/maps/scripts/Update_Map_HMI.py` can be used/referenced to update maps as new data comes in.

### Custom Use Cases

MagMAP is a python library at its core, and while the above two executable scripts are suitable for most general use cases, the library contains additional features and customizations when imported and used in a custom script.  See the `magmap/scripts` folder for examples.

--------------------------------
  
## Contact

James Turtle ([jturtle@predsci.com](mailto:jturtle@predsci.com))

## Contributors

- James Turtle
- Ronald M. Caplan
- Jon A. Linker
- Cooper Downs
