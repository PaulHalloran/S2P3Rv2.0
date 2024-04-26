#!/usr/bin/env python
import time
import process_cmip6_for_s2p3_rv20_improved_interpolation_multiprocessor as proc_cmip
import os
import tempfile

min_depth_lim = 4.0
max_depth_lim = 50.0

start_year = 2015
end_year = 2015



cmip_models = ['UKESM1-0-LL'] # note that this should match exactly the model name used in the filename
# experiments = ['historical','ssp119','ssp585'] # note that this shoudl match exactly the experiment name used in the filename
#experiments = ['historical','ssp585'] # note that this shoudl match exactly the experiment name used in the filename
experiments = ['ssp245'] # note that this shoudl match exactly the experiment name used in the filename
my_suffix = '_r1i1p1f2_all.nc' #e.g. '_all.nc'
my_suffix_windspeed_output = '_r1i1p1f2_all.nc' #e.g. '_all.nc'

# location of source code
base_directory = os.path.join(os.environ["HOME"],"code", "S2P3Rv2.0")

domain_file = 's12_m2_s2_n2_h_map_SundaShelf.dat'
# Note, the script will fail with the error 'OverflowError: cannot serialize a string larger than 2 GiB'
# if this file is too big. For example, a global 4km dataset is too big, but 4km from 30S to 30N is OK
# The limitation is that the data has to be pickled to be run in a parallellised way, and currently
# pickle has a 2GB limit in python2


#Specify where the CMIP data is stored on your computer
base_directory_containing_files_to_process = '/project/ciid/projects/WISERAP_Sunda_Shelf/'
base_output_directory = '/project/ciid/projects/WISERAP_Sunda_Shelf/'

# note Exeter doing temporary stuff on RAMdisk to speed things up
# using TMPDIR here

#local_temp = os.getenv('LOCALTEMP')
#with tempfile.TemporaryDirectory(dir=local_temp) as tmp:

# if using TMPDIR on Slurm don't need to clean up - done by system
base_tmp_output_directory = os.getenv('TMPDIR')+'/'
print(f"temp dir is {base_tmp_output_directory}")

directory_containing_land_sea_mask_files = '/project/champ/data/CMIP6/CMIP/MOHC/' + \
	'UKESM1-0-LL/piControl/r1i1p1f2/fx/sftlf/gn/v20190705/'

# call main programme
start = time.time()
proc_cmip.main(min_depth_lim, max_depth_lim,start_year,end_year,cmip_models,experiments,my_suffix,
	my_suffix_windspeed_output,base_directory, domain_file,base_directory_containing_files_to_process,
	base_output_directory,base_tmp_output_directory,directory_containing_land_sea_mask_files)

end = time.time()
diff = end - start
print(f"elapsed time: {diff}")