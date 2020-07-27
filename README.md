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
```

- CDO (https://code.mpimet.mpg.de/projects/cdo/)
	- instructions: https://code.mpimet.mpg.de/projects/cdo/wiki#Download-Compile-Install note this can be tricky - make sure that it used the same netcdf libraries and hdf libraries as your python bits and pieces

- The OSU Tidal Data Inversion software, installation described below (http://volkov.oce.orst.edu/tides/)

## Installation


clone this repository to your computer

```
git clone https://github.com/PaulHalloran/S2P3Rv2.0.git
```

Move into the S2P3Rv2.0 directory

```
cd S2P3Rv2.0
```

# Setting up and generating forcing data

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

if required install gfortran

sudo apt-get install gfortran

compile the code (using the three lines below), not worrying about 'warning' messages:

gfortran -o extract_HC -fconvert=swap -frecord-marker=4 extract_HC.f90 subs.f90

gfortran -o predict_tide -fconvert=swap -frecord-marker=4 predict_tide.f90 subs.f90

gfortran -o extract_local_model -fconvert=swap -frecord-marker=4 extract_local_model.f90 subs.f90


##################################################################
#  Producing the domain file containing the latitude, longitude, tidal and bathymetry data  #
##################################################################

edit the first 6 lines in tides.py to specify the domain you want the model to run for and the horizonal resolution
(lines copied below for illustration)

minimum_latitude = 5.0
maximum_latitude = 10.0
latitude_resolution = 0.5 #degrees

minimum_longitude = -100.0
maximum_longitude = -90.0
longitude_resolution = 0.5 #degrees

produce the domain file by running tides.py

python2.7 tides.py

You can check the output by opening the file s12_m2_s2_n2_h_map.dat

Note that tides_bathymetry.py is an alternative script for generating the domain file. This gives you the opportunity to specify your own bathymetry file to be used for the bathymetry to the s2p3_r2.0 model rather than using that which is integrated into OTPS2  

##################################################################
#  Producing the meteorological files  option 1, from NCEP       #
##################################################################

*NOTE THESE NCEP INSTRUCTIONS DO NOT INDLUDE DOWNWELING LOMG AND SHORTWAVE USE THE DEPRECIATED MODEL VERSION s2p3_rv2.0_no_prescribed_radiation.f90*

make a directory to hold the meteorological forcing data:

mkdir met_data

The forcing data can be downloaded from NCEP (https://www.esrl.noaa.gov/psd/data/gridded/data.ncep.reanalysis.surface.html). Note that the sig995 variables actually correspond to an atmospheric level which is 99.5% the pressure of that at the surface, whcih typically corresponds to a height of about 42m above ground level. Therefore treat the NCEP forcing example here as indicative of how to run teh model, rather than teh best dataset to use.

Download these data with:

total relative humidity

wget ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.dailyavgs/surface/rhum.sig995.*.nc -P met_data

total sea level pressure

wget ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.dailyavgs/surface/pres.sfc.*.nc -P met_data

2m air temperature

wget ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.dailyavgs/surface_gauss/air.2m.gauss.*.nc -P met_data

total cloud cover

wget ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.dailyavgs/other_gauss/tcdc.eatm.gauss.*.nc -P met_data

wind v

wget ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.dailyavgs/surface/vwnd.sig995.*.nc -P met_data

wind u

wget ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.dailyavgs/surface/uwnd.sig995.*.nc -P met_data

edit process_ncep_for_s2p3_rv2.0.py to specify the spatial resolution you want the atm. forcing data to be at, and the years for which you want to perform the run

produce the final forcing files with:

python2.7 process_ncep_for_s2p3_rv2.0.py

##################################################################
#  Producing the meteorological files  option 2, from CMIP       #
##################################################################

make a directory to hold the meterological forcing data:

mkdir met_data

The forcing data can be downloaded from the cmip5 archive (https://esgf-node.llnl.gov/projects/esgf-llnl/)

The required variables are, U and V surface winds, cloud fraction, specific humidity (from which we calculate relative humidity), surface air temperature, sea-level pressure, net downwelling shortwave, downwelling longwave and land fraction

The CMIP5 variable names for these are:
vas, uas, clt, huss, tas, psl, rsds, rlds, sftlf

Note, the land fraction is used to replace values from atm. grid cells over land with the value from the nearest neighbouring over-ocean cell. This has been implemented to avoid (e.g.) anomalously low wind speeds arising from high terrestrial surface roughness occuring over the sea.

They must all be downloaded at daily frequency. At present the code has been set up only to work with a single ensemble member.

once downloaded, the multiple files for each variable within a model must be merged into a single file with a name in the format MODELNAME_VARIABLENAME_EXPERIMENT_NAME_ENSEMBLENAME.nc. This can be done with cdo, e.g.

cdo mergetime tas*MIROC-ESM_historical_r1i1p1*.nc MIROC-ESM_tas_historical_all.nc

edit process_cmip5_for_s2p3_rv2.0.py to specify the:
 - spatial resolution you want the atm. forcing data to be at
 - the years for which you want to perform the run
 - the name of teh cmip model you want to process (cmip_model = )
 - the experiment name you want to process (experiment = )
 - the location of the merged netcdf files for that model/experiment (directory_containing_files_to_process = )

run the script to produce the forcing with:

python2.7 process_cmip5_for_s2p3_rv2.0.py

##################################################################
#  Producing the meteorological files  option 3, from ECMWF      #
##################################################################

- Create an account with ECMWF and log in
- Go to https://confluence.ecmwf.int/display/WEBAPI/Access+ECMWF+Public+Datasets
- Follow the instructions to set up python script based downloads

For ERA 20C create python scripts like:

This is for net downwards solar radiation at surface, downwards longwave radiation at surface, sea level pressure, 2m air temperature, total cloud cover and 10m U and V winds and relative humidity.

** Note before doing this for shortwave and longwave (parameters 176.128 and 175.128) look at the note at the bottom of this section, this might require some minor changes - suggestion at bottom of section*


```
from ecmwfapi import ECMWFDataServer

server = ECMWFDataServer(url="https://api.ecmwf.int/v1",key="eaf46da7f44cde27c971a213957fd86d",email="p.halloran@exeter.ac.uk")

server.retrieve({
    "class": "e2",
    "dataset": "era20c",
    "date": "1900-01-01/to/2010-12-31",
    "expver": "1",
    "levtype": "sfc",
    "param": "175.128",
    "step": "3",
    "stream": "oper",
    "time": "00:00:00/06:00:00/12:00:00/18:00:00",
    "type": "fc",
    "target": "downwelling_longwave.grb",
		"resol": "auto",
		"grid": "1.0/1.0",
})

server.retrieve({
    "class": "e2",
    "dataset": "era20c",
    "date": "1900-01-01/to/2010-12-31",
    "expver": "1",
    "levtype": "sfc",
    "param": "176.128",
    "step": "3",
    "stream": "oper",
    "time": "00:00:00/06:00:00/12:00:00/18:00:00",
    "type": "fc",
    "target": "downwelling_shortwave.grb",
		"resol": "auto",
		"grid": "1.0/1.0",
})

server.retrieve({
"class": "e2",
"dataset": "era20c",
"date": "1900-01-01/to/2010-12-31",
"expver": "1",
"levtype": "sfc",
"param": "167.128",
"stream": "oper",
"time": "00:00:00/06:00:00/12:00:00/18:00:00",
"type": "an",
"target": "tas.grb",
"resol": "auto",
"grid": "1.0/1.0",
})

server.retrieve({
"class": "e2",
"dataset": "era20c",
"date": "1900-01-01/to/2010-12-31",
"expver": "1",
"levtype": "sfc",
"param": "165.128",
"stream": "oper",
"time": "00:00:00/06:00:00/12:00:00/18:00:00",
"type": "an",
"target": "uas.grb",
"resol": "auto",
"grid": "1.0/1.0",
})

server.retrieve({
"class": "e2",
"dataset": "era20c",
"date": "1900-01-01/to/2010-12-31",
"expver": "1",
"levtype": "sfc",
"param": "166.128",
"stream": "oper",
"time": "00:00:00/06:00:00/12:00:00/18:00:00",
"type": "an",
"target": "vas.grb",
"resol": "auto",
"grid": "1.0/1.0",
})

server.retrieve({
"class": "e2",
"dataset": "era20c",
"date": "1900-01-01/to/2010-12-31",
"expver": "1",
"levtype": "sfc",
"param": "164.128",
"stream": "oper",
"time": "00:00:00/06:00:00/12:00:00/18:00:00",
"type": "an",
"target": "clt.grb",
"resol": "auto",
"grid": "1.0/1.0",
})

server.retrieve({
"class": "e2",
"dataset": "era20c",
"date": "1900-01-01/to/2010-12-31",
"expver": "1",
"levtype": "sfc",
"param": "151.128",
"stream": "oper",
"time": "00:00:00/06:00:00/12:00:00/18:00:00",
"type": "an",
"target": "psl.grb",
"resol": "auto",
"grid": "1.0/1.0",
})

server.retrieve({
"class": "e2",
"dataset": "era20c",
"date": "1900-01-01/to/2010-12-31",
"expver": "1",
"levelist": "1000",
"levtype": "pl",
"param": "157.128",
"stream": "oper",
"time": "12:00:00:00/06:00:00/12:00:00/18:00:00",
"type": "an",
"target": "rh.grb",
"resol": "auto",
"grid": "1.0/1.0",
})

server.retrieve({
"class": "e2",
"dataset": "era20c",
"date": "1900-01-01",
"levtype": "sfc",
"param": "172.128",
"stream": "oper",
"time": "00:00:00",
"target": "lsmask_output.grb",
"resol": "auto",
"grid": "1.0/1.0",
})
```


-run this scripts in python2.7

note, alternatively all surface fields can be downloaded in one file and split with e.g.:

grib_copy -w shortName=2t psl_t_tcc_u_v_output.grb tas_extract.grb


- convert grib files to netcdf:

NOTE for bug files use cdo instead e.g.:

cdo -f nc copy rsds.grib day_mean/rsds_all.nc

or for smaller files:

grib_to_netcdf rh.grib -o day_mean/hurs.nc

grib_to_netcdf total_cloud_cover.grib -o day_mean/clt.nc

grib_to_netcdf mean_sea_level_pressure.grib -o day_mean/psl.nc

grib_to_netcdf 10m_v_component_of_wind.grib -o day_mean/vas.nc

grib_to_netcdf 10m_u_component_of_wind.grib -o day_mean/uas.nc

grib_to_netcdf 2m_temperature.grib -o day_mean/tas.nc

** Note before doing this for shortwave and longwave look at the note at the bottom of this section, this might require some thinking*

grib_to_netcdf downwelling_longwave.grib -o day_mean/rlds.nc

grib_to_netcdf downwelling_shortwave.grib -o day_mean/rsds.nc


convert (e.g.) 3 hourly values to a daily mean and extract the region of interest (the later just helps to make processing faster)

cdo daymean -sellonlatbox,-180,180,-30,30 tas.nc region_of_interest/tas.nc

cdo daymean -sellonlatbox,-180,180,-30,30 clt.nc region_of_interest/clt.nc

cdo daymean sellonlatbox,-180,180,-30,30 hurs.nc region_of_interest/hurs.nc

cdo daymean sellonlatbox,-180,180,-30,30 psl.nc region_of_interest/psl.nc

cdo daymean sellonlatbox,-180,180,-30,30 uas.nc region_of_interest/uas.nc

cdo daymean sellonlatbox,-180,180,-30,30 vas.nc region_of_interest/vas.nc

cdo daymean sellonlatbox,-180,180,-30,30 rlds.nc region_of_interest/rlds.nc

cdo daymean sellonlatbox,-180,180,-30,30 rsds.nc region_of_interest/rsds.nc

I have not yet produced a script to convert this to a forcing file, but I suggest making minor changes to process_ecmwf_era5_for_s2p3_rv2.0.py and using that.
in process_ecmwf_era5_for_s2p3_rv2.0.py these two lines:

znew[:,np.where(input_variables2 == 'rsds')[0],:] /= (60.0*60.0) # ECMWF ERA5 data is in J/m2/hour - see https://confluence.ecmwf.int/display/CKB/ERA5+data+documentation#ERA5datadocumentation-Dataformat
znew[:,np.where(input_variables2 == 'rlds')[0],:] /= (60.0*60.0) # ECMWF ERA5 data is in J/m2/hour - see https://conf

convert the units from J/m2/hour to W/m2

I *think* in all other ECMWF products (including ERA20C) the values are J/m2 accumulated since the start of that day, i.e. a better way to get the mean would be to download just the last value from each day (or maybe it is 1st value of the next day) and divide that b the number of seconds up to that point.

##################################################################
#  ECMWF ERA5                                                    #
##################################################################

Note, ERA 5 is downloaded differently and means are different, from here: https://cds.climate.copernicus.eu

but example retrieval script is:

example_ecmwf_era5_retrieval_script.py

and processing script is:

process_ecmwf_era5_for_s2p3_rv2.0.py

##################################################################
#  Producing the nutrient initialisation file                    #
##################################################################

edit the first three lines of

initialisation_nitrate.py

(hopefully slef explanatory)

then run with:

python2.7 initialisation_nitrate.py


##################################################################
# forcing files                                                   #
##################################################################

The files:

s12_m2_s2_n2_h_map.dat, initial_nitrate.dat

and those listed by

ls met_data/*.dat

are those required to run the model these need to be copied to the /domain and meteorology directories where the model has been set up respectively (see readme for running the model)
