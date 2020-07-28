# S2P3-R v2.0 README

The instructions here assume you are working on a linux machine

##  Introduction

This repository relates to the manuscript 'S2P3-R v2.0: computationally efficient modelling of shelf seas on regional to global scales' (Halloran et al., 2020) submitted to GMD

The repository contains three top level directories (forcing,  model and  processing). The code contained within 'Forcing' generates the forcing files required by the model from user selected input datastes and criteria. The code contained within 'Model' pertains to running the model and generating model output. The code contained within 'Processing' gives some example post-processing and plotting scripts to help the user test the model and produce some standard forms of output. These steps correspond to steps in the flow diagram presented in figure 3 of Halloran et al., (2020).

This README takes the user through the initial requirements, then generating the forcing files, running the model and undertaking basic post-processing.  

##  Requirements

- GFortran compiler (https://gcc.gnu.org/wiki/GFortran)

```
sudo apt-get install gfortran
```

- git (https://gist.github.com/derhuerst/1b15ff4652a867391f03)

	- install with:

```
sudo apt-get update
sudo apt-get upgrade
sudo apt-get install git
```

- Python2.7 with additional libraries (see below)
	- installing conda will make this easier https://conda.io/docs/user-guide/install/index.html
		- additional libraries can then be installed with:

```
conda install pandas
conda install numpy
conda install -c conda-forge iris
conda install -c conda-forge iris-grib
conda install matplotlib
conda install -c conda-forge gridfill
conda install -c conda-forge cdsapi 
```

- CDO (https://code.mpimet.mpg.de/projects/cdo/)
	- instructions: https://code.mpimet.mpg.de/projects/cdo/wiki#Download-Compile-Install note this can be tricky - make sure that it used the same netcdf libraries and hdf libraries as your python bits and pieces

- The OSU Tidal Data Inversion software, installation described below (http://volkov.oce.orst.edu/tides/)

## Installation and 1st setup


clone this repository to your computer

```
git clone https://github.com/PaulHalloran/S2P3Rv2.0.git
```

Move into the S2P3Rv2.0 directory

```
cd S2P3Rv2.0
```

### Setting up and generating forcing data

Tidal slope is an important forcing to the model. The forcing scripts generate this using the Oregon State University TPOX tidal data (https://www.tpxo.net/home)

First ownload OSU data and software:

download the TPXO9.1 bin file

```
wget ftp://ftp.oce.orst.edu/dist/tides/Global/tpxo9.tar.gz
```

uncompress this file:

```
tar zxvf tpxo9.tar.gz
```

get the tidal calculation software

```
wget ftp://ftp.oce.orst.edu/dist/tides/OTPS2.tar.Z
```

uncompress this

```
tar zxvf OTPS2.tar.Z
```

move the data you have downlaoded into the current working dorectory (ignoring the  warning about 'cannot move...'):

```
mv OTPS2/* .
```

Move the data out of the subdirectory to your working directory:

```
mv OTPS2/DATA/load_file DATA/
```

Move our own version of the Model_atlas file (specific to our requirements here) to our data directory

```
mv Model_atlas DATA/
```

Compile the OSU tidal model code (using the three lines below), not worrying about 'warning' messages:

```
gfortran -o extract_HC -fconvert=swap -frecord-marker=4 extract_HC.f90 subs.f90
gfortran -o predict_tide -fconvert=swap -frecord-marker=4 predict_tide.f90 subs.f90
gfortran -o extract_local_model -fconvert=swap -frecord-marker=4 extract_local_model.f90 subs.f90
```

### Producing the domain file
This file contains the latitude, longitude, tidal and bathymetry data for the simulation's domain

- Obtain your chosen bathymetry file in netcdf format. The global simulations in Halloran et al., 2020 use the ETOPO1 product (obtained from https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/bedrock/grid_registered/netcdf/)

edit the first 7 lines in tides_bathymetry.py to specify the domain you want the model to run for and the horizonal resolution
(lines copied below for illustration)

```
minimum_latitude = 5.0
maximum_latitude = 10.0
latitude_resolution = 0.5 #degrees

minimum_longitude = -100.0
maximum_longitude = -90.0
longitude_resolution = 0.5 #degrees

location_of_bathymetry_file =  '/pathway_to_where_you_have_downloaded_your_bathymetry_file/ETOPO1_Bed_g_gmt4.nc'
```

produce the domain file by tides_bathymetry.py

```
python2.7 tides_bathymetry.py
```

You can check the output has been writen correctly by ensuring the file s12_m2_s2_n2_h_map.dat contains non-zero data

###  Producing the meteorological files

The meterological frocing can come from a variety ofatmopsheric models/reanalyses. Here we provide scripts and instructions for generating this forcing from the ECMWF ERA5 product and from CMIP models

####  OPTION 1: Producing the meteorological files from CMIP

make a directory to hold the meterological forcing data:

```
mkdir met_data
```

The forcing data can be downloaded from the cmip5 archive (https://esgf-node.llnl.gov/projects/esgf-llnl/)

The required variables are, U and V surface winds, specific humidity (from which we calculate relative humidity), surface air temperature, sea-level pressure, net downwelling shortwave, downwelling longwave and land fraction

The CMIP5 variable names for these are:
vas, uas, huss, tas, psl, rsds, rlds, sftlf

Note, the land fraction is used to replace values from atm. grid cells over land with the value from the nearest neighbouring over-ocean cell. This has been implemented to avoid (e.g.) anomalously low wind speeds arising from high terrestrial surface roughness occuring over the sea.

They must all be downloaded at daily frequency. At present the code has been set up only to work with a single ensemble member.

once downloaded, the multiple files for each variable within a model must be merged into a single file with a name in the format MODELNAME_VARIABLENAME_EXPERIMENT_NAME_ENSEMBLENAME.nc. This can be done with cdo, e.g.

```
cdo mergetime tas*MIROC-ESM_historical_r1i1p1*.nc MIROC-ESM_tas_historical_all.nc
```

edit process_cmip6_for_s2p3_rv2.0.py or process_cmip5_for_s2p3_rv2.0.py to specify the:
 - spatial resolution you want the atm. forcing data to be at
 - the years for which you want to perform the run
 - the name of the cmip model you want to process (cmip_model = )
 - the experiment name you want to process (experiment = )
 - the location of the merged netcdf files for that model/experiment (directory_containing_files_to_process = )

run the script to produce the forcing with (e.g.):

```
python2.7 process_cmip6_for_s2p3_rv2.0.py
```

####  OPTION 2: Producing the meteorological files from ECMWF's ERA5 reanalysis

ERA5 data is available from here: https://cds.climate.copernicus.eu

Data can be retrieved automatically by running the python script 'example_ecmwf_era5_retrieval_script_netcdf.py':

```
python2.7 example_ecmwf_era5_retrieval_script_netcdf.py
```

Note that you may wish you edit the list of years under each 'year' heading to download more or less data. 

- Once the data is downloaded, iedit the file process_ecmwf_era5_for_s2p3_rv2.0.py to specify the:
      - spatial resolution you want the atm. forcing data to be at
      - the years for which you want to perform the run
      - the name of the cmip model you want to process (cmip_model = )
      - the experiment name you want to process (experiment = )
      - the location of the merged netcdf files for that model/experiment (directory_containing_files_to_process = )

- run 'process_ecmwf_era5_for_s2p3_rv2.0.py'

```
python2.7 process_ecmwf_era5_for_s2p3_rv2.0.py
```

###  Producing the nutrient initialisation file

Download the gridded netcdf version of the World Ocean Atlas nitrate data from here: https://www.nodc.noaa.gov/cgi-bin/OC5/woa13/woa13oxnu.pl
- These instructions assume you are using teh annual mean file: n00_01.nc, but it may be advisable to use the file pertaining to the winter season in the hemisphere of your simulation, because this is what the model is expecting 

Edit the strings on the right hand side of the first three lines of 'initialisation_nitrate.py'. These should point to the chosen output file name for the nutrient file, the name of the domain file produced above, and the location of the World Ocean Atlas data respectively.

output_file_name = 'initial_nitrate.dat'
domain_file_name = 's12_m2_s2_n2_h_map.dat'
location_of_World_Ocean_Atlas13_Nitrate_file =  '/the_path_to_your_downloaded_World_ocean_atlas_data/woa13_all_n00_01.nc'

generate teh nitrate forcing ancillary with:

```
python2.7 initialisation_nitrate.py
```

### Summary of the forcing files

The files:

s12_m2_s2_n2_h_map.dat, initial_nitrate.dat

and those listed by

```
ls met_data/*.dat
```

are those required to run the model these need to be copied to the /domain and meteorology directories where the model has been set up respectively (see readme for running the model)

## Setting up and running the model


*Hint. You may want to speed things up by creating a RAM disk and making this the temporary location to hold the unzipped met data.

e.g.
```
sudo mkdir /mnt/ramdisk
sudo mount -t tmpfs -o size=3g tmpfs /mnt/ramdisk
then your temp location is  /mnt/ramdisk:
```
*

Move into the 'model' directory within S2P3Rv2.0

cd S2P3Rv2.0/model

- compile code

```
gfortran -Ofast -o s2p3_rv2.0 s2p3_rv2.0.f90
```

s2p3_rv2.0 should now be a working executable

- copying in the forcing data

We want to copy the forcing files created above and listed in the section 'Summary of the forcing files' to the the location from which the model will access them. This may be a different machine if running on a compute server or HPC.

If working on a simple computer and following teh instructions exactly as above this would involve:

*Where **my_path** is the path to the directory where you cloned the repository from github, s12_m2_s2_n2_h_map.dat is the domain file you produced in the sectoin titled 'Producing the domain file', initial_nitrate.dat is the file you produced in the section titled 'Producing the nutrient initialisation file' and the files starting with the name 'met_data_' are those produced under the section titled 'Producing the meteorological files'.*

*Where **my_meterology_path** is location you specified as the 'output_directory' in the script with a name like 'process_x_for_s2p3_rv2.0.py', and **my_met_path_for_model_runs** is the location you want the model to read the files from.*

*Note that there is no need to copy these files, you could just point to the location where they were produced if not undertaking muyltiple different runs.*

```
cp /my_path/S2P3Rv2.0/forcing/s12_m2_s2_n2_h_map.dat /my_path/S2P3Rv2.0/model/domain/
cp /my_path/S2P3Rv2.0/forcing/initial_nitrate.dat /my_path/S2P3Rv2.0/model/domain/
cp -r my_meterology_path/met_data_*.tar.gz my_met_path_for_model_runs/
 
```

scp example_location_1/s2p3_rv2.0_forcing/initial_nitrate.dat example_location_2/s2p3_rv2.0/domain

scp -r example_location_1/s2p3_rv2.0_met_data/met_data_*.tar.gz example_location_2/my_runs_met_data

#running the model

before running you will need to edit at least one line in:

s2p3_rv2.0/main/run_map_parallel.py

Here you must change the line 'base_directory = ...' to point to the s2p3_rv2.0 directory on your computer

You will also probably have to edit the 'met_data_location = ' to point to the locatino you have copied the tar.gz met files to (example_location_2/my_runs_met_data in the example above)

You may also have to create the tempoarry directory which will hold the untarred and gziped met files (e.g. met_data_temporary_location = base_directory+'met/spatial_data/')

You may well also want to change the number or processors, output file names etc. here

Under the heading 'Variables to output from model' you can select which variables to output from the model, see the model output section below.

then run with either:

python2.7 run_map_parallel.py

OR if you are running on a cluster/supercomputer with a batch system, you may be able to use a runscript like runscript_parallel and submit with something like:

msub runscript_parallel

#model output

The model output will be in the directory specified wiuthin run_map_parallel.py, with the default output file being output_map
This contains in columns, left to right: day number, longitude, latitude, then the variables you have specified in run_map_parallel.py in the order specified under the heading 'Variables to output from model'.

Note changing this will necessitate changes to the processing into netcdf script here https://bitbucket.org/paulhalloran/s2p3_rv2.0_processing


############################################################
# Processing/plotting the output                           #
############################################################

see https://bitbucket.org/paulhalloran/s2p3_rv2.0_processing

Scripts and instructions are provided to convert the output into netcdf format and plot the output


############################################################
# Points to note                                           #
############################################################

############################################################
# Updates                                                  #
############################################################

To avoid overly hot SSTs in some tropical areas, changes being made to prescribe downwelling shortwave and longwave radiation. The net downward shortwave at the surface is used directly by the model, but the longwave requires a minor change to the model code. I've changed:
hl=0.985*5.67d-8*(surf_temp**4.0)*(0.39-0.05*vap**0.5)*(1.0-0.6d-4*cloud(idmet)**2.0)   !original codes Longwave
to
hl=0.985*5.67d-8*(surf_temp**4.0)-lw_rad_down0  !Longwave
note, I'm not sure why humiditry would have been used in tgeh original calculation. Is this a way at getting towards skin temperature from tghe model's surface level temperature, or accounting for absorbtion by non-cloud water vapour?
