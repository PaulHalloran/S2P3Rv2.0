# S2P3-R v2.0 README

The instructions here assume you are working on a linux machine

##  Introduction

This repository relates to the manuscript 'S2P3-R v2.0: computationally efficient modelling of shelf seas on regional to global scales' (Halloran et al., 2020) submitted to GMD

The repository contains three top level directories (forcing,  model and  processing). The code contained within 'Forcing' generates the forcing files required by the model from user selected input datasets and criteria. The code contained within 'Model' pertains to running the model and generating model output. The code contained within 'Processing' gives some example post-processing and plotting scripts to help the user test the model and produce some standard forms of output. These steps correspond to steps in the flow diagram presented in figure 3 of Halloran et al., (2020).

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
**NOTE from this point in the instructions *my_path* refers to the path specifying the directory where you have cloned the repository from github**

Move into the S2P3Rv2.0 directory

```
cd S2P3Rv2.0
```

### Setting up and generating forcing data

Move into the forcing directory

```
cd forcing
```

Tidal slope is an important forcing to the model. The forcing scripts generate this using the Oregon State University TPOX tidal data (https://www.tpxo.net/home)

To get hold of this you now need to register by following the instructions here: https://www.tpxo.net/tpxo-products-and-registration
The files you require are grid_tpxo9 and u_tpxo9.v1
Available within the University of Exeter at https://universityofexeteruk-my.sharepoint.com/:f:/g/personal/p_halloran_exeter_ac_uk/ElYRNQ1oNllArEiHkYU5bSkByTXvNSieTddNXzuK9oTpGA?email=P.Halloran%40exeter.ac.uk&e=RpeqCP


get the tidal calculation software: https://drive.google.com/file/d/1FBlS_Xmf6_dnCg1T0t5GSTRTwMjLuA8N/view

uncompress this

```
tar zxvf OTPS.tar.Z
```

move the data you have downloaded into the current working directory:

```
mv OTPS/* .
```



Move our own version of the Model_atlas file (specific to our requirements here) to our data directory

```
mv Model_atlas DATA/
```

Copy the files grid_tpxo9 and u_tpxo9.v1 obtained above into the DATA directory

Compile the OSU tidal model code (using the three lines below), not worrying about 'warning' messages:

```
gfortran -o extract_HC -fconvert=swap -frecord-marker=4 extract_HC.f90 subs.f90
gfortran -o predict_tide -fconvert=swap -frecord-marker=4 predict_tide.f90 subs.f90
gfortran -o extract_local_model -fconvert=swap -frecord-marker=4 extract_local_model.f90 subs.f90
```

### Producing the domain file
This file contains the latitude, longitude, tidal and bathymetry data for the simulation's domain

- Obtain your chosen bathymetry file in netcdf format. The global simulations in Halloran et al., 2020 use the ETOPO1 product ETOPO1_Bed_g_gmt4.nc (obtained from https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/bedrock/grid_registered/netcdf/)

edit the first 7 lines in tides_bathymetry.py to specify the domain you want the model to run for and the horizontal resolution
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

You can check the output has been written correctly by ensuring the file s12_m2_s2_n2_h_map.dat contains non-zero data

###  Producing the meteorological files

The meteorological forcing can come from a variety of atmospheric models/reanalyses. Here we provide scripts and instructions for generating this forcing from the ECMWF ERA5 product and from CMIP models

####  OPTION 1: Producing the meteorological files from CMIP

make a directory within S2P3Rv2.0 to hold the meteorological forcing data:

```
mkdir ../met_data
```

The forcing data can be downloaded from the CMIP6 archive (https://esgf-node.llnl.gov/projects/esgf-llnl/)

The required variables are, U and V surface winds, specific humidity (from which we calculate relative humidity), surface air temperature, sea-level pressure, net downwelling shortwave, downwelling longwave and land fraction

The CMIP variable names for these are:
vas, uas, hurs, tas, psl, rsds and rlds at daily frequency and sftlf as a fixed frequency file

Note, the land fraction is used to replace values from atmospheric grid cells over land with the value from the nearest neighbouring over-ocean cell. This has been implemented to avoid (e.g.) anomalously low wind speeds arising from high terrestrial surface roughness occurring over the sea.

They must all be downloaded at daily frequency. At present the code has been set up only to work with a single ensemble member.

once downloaded, the multiple files for each variable within a model must be merged into a single file with a name in the format MODELNAME_VARIABLENAME_EXPERIMENT_NAME_ENSEMBLENAME.nc. This can be done with cdo, e.g.

```
cdo mergetime tas*MIROC-ESM_historical_r1i1p1*.nc tas_MIROC-ESM_historical_all.nc
```

edit process_cmip6_for_s2p3_rv2.0.py to specify the:
 - spatial resolution you want the atm. forcing data to be at
 - the years for which you want to perform the run
 - the name of the CMIP model you want to process (cmip_model = )
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

- Once the data is downloaded, edit the file process_ecmwf_era5_for_s2p3_rv2.0.py to specify the:
      - spatial resolution you want the atm. forcing data to be at
      - the years for which you want to perform the run
      - the name of the CMIP model you want to process (cmip_model = )
      - the experiment name you want to process (experiment = )
      - the location of the merged netcdf files for that model/experiment (directory_containing_files_to_process = )

- run 'process_ecmwf_era5_for_s2p3_rv2.0.py'

```
python2.7 process_ecmwf_era5_for_s2p3_rv2.0.py
```

###  Producing the nutrient initialisation file

Download the gridded netcdf version of the World Ocean Atlas nitrate data from here: https://www.nodc.noaa.gov/cgi-bin/OC5/woa13/woa13oxnu.pl
- These instructions assume you are using the annual mean file: n00_01.nc, but it may be advisable to use the file pertaining to the winter season in the hemisphere of your simulation, because this is what the model is expecting

Edit the strings on the right hand side of the first three lines of 'initialisation_nitrate.py'. These should point to the chosen output file name for the nutrient file, the name of the domain file produced above, and the location of the World Ocean Atlas data respectively.

output_file_name = 'initial_nitrate.dat'
domain_file_name = 's12_m2_s2_n2_h_map.dat'
location_of_World_Ocean_Atlas13_Nitrate_file =  '/the_path_to_your_downloaded_World_ocean_atlas_data/woa13_all_n00_01.nc'

Generate the nitrate forcing ancillary with:

```
python2.7 initialisation_nitrate.py
```

### Summary of the forcing files

The files:

s12_m2_s2_n2_h_map.dat, initial_nitrate.dat

and unless you have specified a different output directory for the meteorological data those listed by:

```
ls met_data/*.tar.gz
```

are those required to run the model these need to be copied to the /domain and meteorology directories where the model has been set up respectively (see readme for running the model)

Note that the meteorological data for each lat/lon location for a specific year exists as a .dat file which is compressed into a single tar.gz file for each year. These are extracted when the model runs, but can be extracted manually to assess their contents.

*For reference the columns in the meteorological files are: day number in year, wind speed (m/s), wind direction (degrees), surface atmospheric temperature (tas, deg C),surface atmospheric temperature duplicated for legacy reasons, sea level pressure (psl, hPa),surface level relative humidity (hurs, %), shortwave downwards radiation at the surface (rsds, wm-2), longwave downwards radiation at the surface (rlds, wm-2).

## Setting up and running the model



*Hint. You may want to speed things up by creating a RAM disk and making this the temporary location to hold the unzipped met data.*

e.g.
```
sudo mkdir /mnt/ramdisk
sudo mount -t tmpfs -o size=3g tmpfs /mnt/ramdisk
```
then your temp location is  /mnt/ramdisk

Move into the 'model' directory within S2P3Rv2.0

```
cd my_path/S2P3Rv2.0/model/main
```

- compile code

```
gfortran -Ofast -o s2p3_rv2.0 s2p3_rv2.0.f90
```

s2p3_rv2.0 should now be a working executable

- copying in the forcing data

We want to copy the forcing files created above and listed in the section 'Summary of the forcing files' to the the location from which the model will access them. This may be a different machine if running on a compute server or HPC.

If working on a simple computer and following the instructions exactly as above this would involve:

*Where **my_path** is the path to the directory where you cloned the repository from github, s12_m2_s2_n2_h_map.dat is the domain file you produced in the sectoin titled 'Producing the domain file', initial_nitrate.dat is the file you produced in the section titled 'Producing the nutrient initialisation file' and the files starting with the name 'met_data_' are those produced under the section titled 'Producing the meteorological files'.*

*Where **my_meterology_path** is location you specified as the 'output_directory' in the script with a name like 'process_x_for_s2p3_rv2.0.py', and **my_met_path_for_model_runs** is the location you want the model to read the files from.*

*Note that there is no need to copy these files, you could just point to the location where they were produced if not undertaking multiple different runs.*

```
cp /my_path/S2P3Rv2.0/forcing/s12_m2_s2_n2_h_map.dat /my_path/S2P3Rv2.0/model/domain/
cp /my_path/S2P3Rv2.0/forcing/initial_nitrate.dat /my_path/S2P3Rv2.0/model/domain/
cp -r my_meterology_path/met_data_*.tar.gz my_met_path_for_model_runs/

```

#running the model

before running you will need to edit the top section of:

```
/my_path/S2P3Rv2.0/model/main/run_map_parallel.py
```

The lines to edit are

```
base_directory = '/my_path/s2p3_rv2.0/' # This is the directory containing the 'forcing', 'model' and 'met' directories.
num_procs = mp.cpu_count() # this will use all available processors. Note that on a multi-node machine the model can only use the processors on one node
# num_procs = 1 # The default is to use all available processors, but it is possible to specify a smaller number of processors.

output_directory = '/some_directory/'  # the directory where you want the output from the model to be written

output_file_name = 'a_filename_to_identify_the_output_from_this_specific_simulation'
meterological_file_name = 'meterological_data' # leave this as it is unless you have changed the code described within the 'Producing the meteorological files' section
domain_file_name = 's12_m2_s2_n2_h_map.dat' # This is the name of the output fine produced by running 'tides_bathymetry.py'
nutrient_file_name = 'initial_nitrate.dat' # This is the name of the output fine produced by running 'initialisation_nitrate.py'
executable_file_name = 's2p3_rv2.0' # The is the compiled model executable, i.e. teh righthand side of the line 'gfortran -Ofast -o s2p3_rv2.0 s2p3_rv2.0.f90' run above.

met_data_location = '/my_met_path_for_model_runs/' # The location containing the tar.gz met files (in the format met_data_year.tar.gz). See the line 'cp -r my_meterology_path/met_data_*.tar.gz my_met_path_for_model_runs/' above

met_data_temporary_location = '/mnt/ramdisk/' # The model uncompresses the meteorological data files to a location from which it can be read quickly. The example here is a RAMdisk (see above), but it can be any storage - ideally fast storage.

start_year = 1950 # The year for which to start the model simulation (note, this should fall within the years for which you have created the meteorological data)
end_year = 2100 # The last year of the simulation. It is is the same as start year the model will simulate for 1 full year
depth_min = 4 # The most shallow water depth to run the model in. NOTE that these numbers MUST be the same as those used in the scripts used to produce the meteorology (e.g. process_cmip6_for_s2p3_rv2.0.py) and nutrient files, otherwise data will not be taken for the correct lats/lons and/or the script will fail
depth_max = 50 # The deepest water depth to run the model for (i.e. so you are not running the model in the open ocean)
write_error_output = False # Change to True for debugging

parallel_processing = True # True if you want to run on more than one processor. A single processor may make some debugging easier.

generate_netcdf_files = True #If True, saves model output as netcdf files. Set to False if you have set write_error_output to True
```

You will also find a list of variable names under the heading 'Variables to output from model' in run_map_parallel.py. Set these to 1 if you want to output this variable, or 0 if you do not wish to output this variable.

- Making sure that you are in the '/my_path/s2p3_rv2.0/model/main/' directory run the model with either:

```
python2.7 run_map_parallel.py
```

OR if you are running on a cluster/supercomputer you may need to submit this with a runscript specific to your batch system. An example using msub is provided in the file 'runscript_parallel'. This would be submitted with 'msub runscript_parallel'

#model output

The model output will be in the directory specified for the 'output_directory' variable in 'run_map_parallel.py'.

If you have specified netcdf output you will have a file for each year and each specified variable.

If you have chosen not to output netcdf files the model will generate a csv file for each year containing form left to right columns of day number, longitude, latitude, then the variables you have specified in run_map_parallel.py in the order specified under the heading 'Variables to output from model'.


## plotting the output
These instructions assume that you have produced output as netcdf files

A single script ('processing/basic_plots.py') is supplied which provides an example of how to read the model output into python, and multiple examples of how to process and plot the data. These examples include:
- Plot a contourmap of all data averaged along the time dimension
- Extract data falling between certain years
- Extract data falling between certain latitude/longitude bounds
- Create a timeseries from mapped data by performing a weighted area average, then plot
- Concert daily data into monthly, seasonal and annual averaged data

At the very least edit the three lines below the line 'Edit the three lines...' to point the script to the model output you wish to plot, then run with:
```
python2.7 -i /my_path/s2p3_rv2.0/processing/vim basic_plots.py
```
