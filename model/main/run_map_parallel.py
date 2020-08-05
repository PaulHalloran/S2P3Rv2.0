
import os
import subprocess
import shutil
import glob
from math import cos, asin, sqrt
import multiprocessing as mp
from functools import partial
import uuid
import time

##################################################
# you may need to change things here             #
##################################################

base_directory = '/home/ph290/s2p3Rv2_manuscript/S2P3Rv2.0/'
num_procs = mp.cpu_count() # this will use all available processors. Note that on a multi-node machine it can only use the processors on one node
# num_procs = 1 # The default is to use all available processors, but it is possible to specify the number of processors.

output_directory = '/data/local_ssd/ph290/s2p3_rv2.0/output/test/'  #where you want the output to go

output_file_name = 'test_output'
meterological_file_name = 'meterological_data'
domain_file_name = 's12_m2_s2_n2_h_map.dat'
nutrient_file_name = 'initial_nitrate.dat'
executable_file_name = 's2p3_rv2.0'

met_data_location = '/data/local_ssd/ph290/test/cmip6/processed/' # The location containing the tar.gz met files (in the format met_data_year.tar.gz)

met_data_temporary_location = '/mnt/ramdisk/' # The location that met data for each year will be un tar.gziped into
# each grid point each year has to read in a new meterology dataset from disk so it may make sense to make this temporary location a RAM disk (see readme)

start_year = 2020

end_year = 2021 # same as start year resuls in a 1 year run
depth_min = 10.0 # NOTE that these numbers MUST be the same as those used in the scripts used to produce the meterology and nutrient files, otherwse data will not be taken for teh correct lats/lons and/or the script will fail
depth_max = 30.0
write_error_output = False

parallel_processing = True

generate_netcdf_files = True
#note does not output error data if write_error_output set to True

#######################################################
# Variables to output from model                      #
# =1 means output, =0 means do not output             #
# the columns in the output will be left to right in  #
# the order that the variables are here top-to-bottom #
#######################################################

include_depth_output=1
include_temp_surface_output=1
include_temp_bottom_output=1
include_chlorophyll_surface_output=1
include_phyto_biomass_surface_output=0
include_phyto_biomass_bottom_output=0
include_PAR_surface_output=1 # daily mean surface of PAR W m-2
include_PAR_bottom_output=1 # daily mean bottom of PAR W m-2
include_windspeed_output=0 # windspeed
include_stressx_output=0 # x component of surface wind drag
include_stressy_output=0 # y component of surface wind drag
include_Etide_output=0 # Mixing power in the tidal currents (assumes constant mixing efficiency 0.003)
include_Ewind_output=0 # Mixing power in the wind (assumes constant mixing efficiency 0.023 and slippage factor=0.025)
include_u_mean_surface_output=0 # daily mean surface u in cm/s !!
include_u_mean_bottom_output=0 # daily mean bottom u in cm/s !!
include_grow1_mean_surface_output=0 # daily mean growth rate d-1 surface
include_grow1_mean_bottom_output=0 # daily mean growth rate d-1 bottom
include_uptake1_mean_surface_output=0 # daily mean DIN uptake rate mmol DIN (mg C)-1 d-1 surface
include_uptake1_mean_bottom_output=0 # daily mean DIN uptake rate mmol DIN (mg C)-1 d-1 bottom
include_tpn1_output=0 # total water column net production / mg C m-2 d-1
include_tpg1_output=0 # total water column gross production / mg C m-2 hd-1
include_speed3_output=0	# depth-mean current speed

columns = [include_depth_output,include_temp_surface_output,include_temp_bottom_output,include_chlorophyll_surface_output,include_phyto_biomass_surface_output,include_phyto_biomass_bottom_output,include_PAR_surface_output,include_PAR_bottom_output,include_windspeed_output,include_stressx_output,include_stressy_output,include_Etide_output,include_Ewind_output,include_u_mean_surface_output,include_u_mean_bottom_output,include_grow1_mean_surface_output,include_grow1_mean_bottom_output,include_uptake1_mean_surface_output,include_uptake1_mean_bottom_output,include_tpn1_output,include_tpg1_output,include_speed3_output]
column_names_all = ['depth','surface temperature','bottom temperature','surface chlorophyll','surface phyto biomass','bottom phyto biomass','surface PAR','bottom PAR','windspeed','stress_x','stress_y','Etide','Ewind','u_mean_surface','u_mean_bottom','grow1_mean_surface','grow1_mean_bottom','uptake1_mean_surface','uptake1_mean_bottom','tpn1','tpg1','speed3']



##################################################
# functions used by the script                   #
##################################################

if generate_netcdf_files:
    import numpy as np
    import iris
    import pandas as pd
    from itertools import compress
    from cf_units import Unit
    column_names = ['day','longitude','latitude']+list(compress(column_names_all, map(bool,columns)))
    specifying_names = False
    ## If specifying_names above is set to True, specify the below. If not, ignore ##
    standard_name=['sea_surface_temperature','sea_surface_temperature','sea_surface_temperature','sea_surface_temperature','sea_surface_temperature']
    long_name=['Sea Surface Temperature','Sea Surface Temperature','Sea Surface Temperature','Sea Surface Temperature','Sea Surface Temperature']
    var_name=['tos','tos','tos','tos','tos']
    units=['K','K','K','K','K']
    if not(specifying_names):
        standard_name=np.tile('sea_surface_temperature',len(column_names))
        long_name=np.tile('Sea Surface Temperature',len(column_names))
        var_name=np.tile('tos',len(column_names))
        units=np.tile('K',len(column_names))



def put_data_into_cube(df,df_domain,variable,specifying_names,standard_name,long_name,var_name,units,run_start_date):
    latitudes = np.unique(df_domain['lat'].values)
    longitudes = np.unique(df_domain['lon'].values)
    latitudes_run = np.unique(df['latitude'].values)
    longitudes_run = np.unique(df['longitude'].values)
    times = np.unique(df['day'].values)
    latitude = iris.coords.DimCoord(latitudes, standard_name='latitude', units='degrees')
    longitude = iris.coords.DimCoord(longitudes, standard_name='longitude', units='degrees')
    # time = iris.coords.DimCoord(times, standard_name='time', units='days')
    time = iris.coords.DimCoord(times, standard_name='time', units=Unit('days since '+run_start_date+' 00:00:0.0', calendar='gregorian'))
    if specifying_names:
        cube = iris.cube.Cube(np.full((times.size,latitudes.size, longitudes.size),-999.99, np.float32),standard_name=standard_name, long_name=long_name, var_name=var_name, units=units,dim_coords_and_dims=[(time,0), (latitude, 1), (longitude, 2)])
    else:
        cube = iris.cube.Cube(np.full((times.size,latitudes.size, longitudes.size),-999.99, np.float32),standard_name=None, long_name=None, var_name=None, units=None,dim_coords_and_dims=[(time,0), (latitude, 1), (longitude, 2)])
    # Z,X,Y = np.meshgrid(cube.coord('time').points,cube.coord('longitude').points,cube.coord('latitude').points)
    data = cube.data.copy()
    # data[:] = -999.99
    days = np.unique(df['day'].values)
    # shape = [X.shape[0],X.shape[2]]
    for i,day in enumerate(days):
        df_tmp = df.loc[df['day'] == day]
        for j,lat in enumerate(df_tmp['latitude'].values):
            lon = df_tmp['longitude'].values[j]
            lat_loc = np.where(np.around(cube.coord('latitude').points,decimals=6) == np.around(lat,decimals=6))[0][0]
            lon_loc = np.where(np.around(cube.coord('longitude').points,decimals=6) == np.around(lon,decimals=6))[0][0]
            data[i,lat_loc,lon_loc] = df_tmp[variable].values[j]
    data = np.ma.masked_where((data.data < -999.9) & (data.data > -1000.0),data)
    # data = np.ma.masked_where(data.data == 0.0,data)
    cube.data = data
    cube.data.fill_value = -999.99
    cube.data.data[~(np.isfinite(cube.data.data))]=-999.99
    cube.data = np.ma.masked_where(cube.data == -999.99,cube.data)
    return cube

def output_netcdf(year,column_names,df,df_domain,specifying_names,standard_name,long_name,var_name,units,run_start_date, output_directory,output_file_name,i):
    column_name = column_names[i]
    output_cube = put_data_into_cube(df,df_domain,column_name,specifying_names,standard_name,long_name,var_name,units,run_start_date)
    iris.fileformats.netcdf.save(output_cube, output_directory+output_file_name+'_'+column_name.replace(" ", "")+'_'+str(year)+'.nc', zlib=True, complevel=2)
    return output_directory+output_file_name+'_'+column_name.replace(" ", "")+'_'+str(year)+'.nc written'


domain_file_for_run = base_directory+'domain/'+domain_file_name
fwidths=[8,8,6,6,6,6,6,6,6,6,6,6,8]
df_domain = pd.read_fwf(domain_file_for_run,names=['lon','lat','t1','t2','t3','t4','t5','t6','t7','t8','t9','t10','depth'],widths = fwidths,
                 skiprows=[0],dtype={'lon':float,'lat':float,'t1':float,'t2':float,'t3':float,'t4':float,'t5':float,'t6':float,'t7':float,'t8':float,'t9':float,'t10':float,'depth':float},usecols=['lon','lat','depth'])


f=open(base_directory+'domain/'+domain_file_name)
lines=f.readlines()
f2=open(base_directory+'domain/'+nutrient_file_name)
lines2=f2.readlines()
lat_domain=[]
lon_domain=[]
alldepth=[]
smaj1=[]
smin1=[]
smaj2=[]
smin2=[]
smaj3=[]
smin3=[]
smaj4=[]
smin4=[]
smaj5=[]
smin5=[]
woa_nutrient=[]
counter = 0
for i,line in enumerate(lines[1::]):
    depth = float(line[77:84])
    if ((depth >= depth_min) & (depth <= depth_max) & (depth > 0.0)):
        lon_domain.append(line[0:8])
        lat_domain.append(line[8:16])
        alldepth.append(line[77:84])
        smaj1.append(line[16:22])
        smin1.append(line[22:28])
        smaj2.append(line[28:34])
        smin2.append(line[34:40])
        smaj3.append(line[40:46])
        smin3.append(line[46:52])
        smaj4.append(line[52:58])
        smin4.append(line[58:64])
        smaj5.append(line[64:70])
        smin5.append(line[70:76])
        woa_nutrient.append(lines2[counter][16:22])
	counter += counter





def run_model(domain_file_name,lats_lons,year,start_year,unique_job_id,met_data_temporary_location,lon_domain,lat_domain,smaj1,smin1,smaj2,smin2,smaj3,smin3,smaj4,smin4,smaj5,smin5,woa_nutrient,alldepth,include_depth_output,include_temp_surface_output,include_temp_bottom_output,include_chlorophyll_surface_output,include_phyto_biomass_surface_output,include_phyto_biomass_bottom_output,include_PAR_surface_output,include_PAR_bottom_output,include_windspeed_output,include_stressx_output,include_stressy_output,include_Etide_output,include_Ewind_output,include_u_mean_surface_output,include_u_mean_bottom_output,include_grow1_mean_surface_output,include_grow1_mean_bottom_output,include_uptake1_mean_surface_output,include_uptake1_mean_bottom_output,include_tpn1_output,include_tpg1_output,include_speed3_output,i):
    #modifying so that the fortran code looks for the correct met file, rather than us having to copy it into the working directory
    # lon,lat = return_domain_lon(base_directory+'domain/'+domain_file_name,i)
    lon_domain_tmp = float(lon_domain[i])
    if lon_domain_tmp < 0.0:
        lon_domain_tmp = 360.0+lon_domain_tmp
    run_command = '\n'.join(['./{} << EOF'.format(executable_file_name),
    str(start_year),
    str(year),
    str(float(lat_domain[i])),
    str(lon_domain_tmp),
    '../domain/{}'.format(domain_file_name),
    '../domain/{}'.format(nutrient_file_name),
    unique_job_id,
    met_data_temporary_location,
    'map',
    str(i+1),
    str(smaj1[i]),
    str(smin1[i]),
    str(smaj2[i]),
    str(smin2[i]),
    str(smaj3[i]),
    str(smin3[i]),
    str(smaj4[i]),
    str(smin4[i]),
    str(smaj5[i]),
    str(smin5[i]),
    str(woa_nutrient[i]),
    str(alldepth[i]),
    str(include_depth_output),
    str(include_temp_surface_output),
    str(include_temp_bottom_output),
    str(include_chlorophyll_surface_output),
    str(include_phyto_biomass_surface_output),
    str(include_phyto_biomass_bottom_output),
    str(include_PAR_surface_output),
    str(include_PAR_bottom_output),
    str(include_windspeed_output),
    str(include_stressx_output),
    str(include_stressy_output),
    str(include_Etide_output),
    str(include_Ewind_output),
    str(include_u_mean_surface_output),
    str(include_u_mean_bottom_output),
    str(include_grow1_mean_surface_output),
    str(include_grow1_mean_bottom_output),
    str(include_uptake1_mean_surface_output),
    str(include_uptake1_mean_bottom_output),
    str(include_tpn1_output),
    str(include_tpg1_output),
    str(include_speed3_output),
    str(start_year),
    str(year),
    str(float(lat_domain[i])),
    str(lon_domain_tmp),
    '../domain/{}'.format(domain_file_name),
    '../domain/{}'.format(nutrient_file_name),
    unique_job_id,
    met_data_temporary_location,
    'map',
    str(i+1),
    str(smaj1[i]),
    str(smin1[i]),
    str(smaj2[i]),
    str(smin2[i]),
    str(smaj3[i]),
    str(smin3[i]),
    str(smaj4[i]),
    str(smin4[i]),
    str(smaj5[i]),
    str(smin5[i]),
    str(woa_nutrient[i]),
    str(alldepth[i]),
    str(include_depth_output),
    str(include_temp_surface_output),
    str(include_temp_bottom_output),
    str(include_chlorophyll_surface_output),
    str(include_phyto_biomass_surface_output),
    str(include_phyto_biomass_bottom_output),
    str(include_PAR_surface_output),
    str(include_PAR_bottom_output),
    str(include_windspeed_output),
    str(include_stressx_output),
    str(include_stressy_output),
    str(include_Etide_output),
    str(include_Ewind_output),
    str(include_u_mean_surface_output),
    str(include_u_mean_bottom_output),
    str(include_grow1_mean_surface_output),
    str(include_grow1_mean_bottom_output),
    str(include_uptake1_mean_surface_output),
    str(include_uptake1_mean_bottom_output),
    str(include_tpn1_output),
    str(include_tpg1_output),
    str(include_speed3_output),
    'EOF'
    ])
    # print run_command
    proc = subprocess.Popen([run_command],  shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = proc.communicate()
    # return out
    return out,err


##################################################
# main program                                   #
##################################################

unique_job_id = str(uuid.uuid4())

num_lines = sum(1 for line in open(base_directory+'domain/'+domain_file_name)) - 1
# num_lines = 10

print 'initial unzipping of met data to extract lat and lon data'
subprocess.call('tar -C '+met_data_temporary_location+' -zxf '+met_data_location+'met_data_'+str(start_year)+'.tar.gz', shell=True)

files = glob.glob(met_data_temporary_location+'*_'+str(start_year)+'.dat')
w, h = 2, len(files) ;
lats_lons = [[0 for x in range(w)] for y in range(h)]
for i,file in enumerate(files):
    tmp = file.split('lat')[-1].split('.dat')[0].split('lon')
    lats_lons[i][0] = float(tmp[0])
    lats_lons[i][1] = float(tmp[1].split('_')[0])



# ##testing
# year=2006
# pool = mp.Pool(processes=num_procs)
# func = partial(run_model, domain_file_name, lats_lons, year, start_year, unique_job_id, met_data_temporary_location,lon_domain,lat_domain,smaj1,smin1,smaj2,smin2,smaj3,smin3,smaj4,smin4,smaj5,smin5,woa_nutrient,alldepth,include_depth_output,include_temp_surface_output,include_temp_bottom_output,include_chlorophyll_surface_output,include_phyto_biomass_surface_output,include_phyto_biomass_bottom_output,include_PAR_surface_output,include_PAR_bottom_output,include_windspeed_output,include_stressx_output,include_stressy_output,include_Etide_output,include_Ewind_output,include_u_mean_surface_output,include_u_mean_bottom_output,include_grow1_mean_surface_output,include_grow1_mean_bottom_output,include_uptake1_mean_surface_output,include_uptake1_mean_bottom_output,include_tpn1_output,include_tpg1_output,include_speed3_output)
# results, errors = zip(*pool.map(func, range(1)))
# print results[0].split('\n')[0]
# print results[0].split('\n')[1]
# print results[0].split('\n')[2]
# # print errors[0].split('\n')[0]
# pool.close()


print 'looping through years'
for year in range(start_year,end_year+1):
# year = start_year
    print year
    #clean up and prexisting met files
    try:
        files_to_delete = glob.glob(met_data_temporary_location+'*.dat')
        [os.remove(f) for f in files_to_delete]
    except:
        print 'no met files to clean up'

    subprocess.call('tar -C '+met_data_temporary_location+' -zxf '+met_data_location+'met_data_'+str(year)+'.tar.gz', shell=True)
    try:
        shutil.move(output_directory+output_file_name+'_'+str(year), output_directory+output_file_name+'_'+str(year)+'_previous')
    except:
        print 'no previous output file to move'

    if parallel_processing:
        pool = mp.Pool(processes=num_procs)
        func = partial(run_model, domain_file_name, lats_lons, year, start_year, unique_job_id, met_data_temporary_location,lon_domain,lat_domain,smaj1,smin1,smaj2,smin2,smaj3,smin3,smaj4,smin4,smaj5,smin5,woa_nutrient,alldepth,include_depth_output,include_temp_surface_output,include_temp_bottom_output,include_chlorophyll_surface_output,include_phyto_biomass_surface_output,include_phyto_biomass_bottom_output,include_PAR_surface_output,include_PAR_bottom_output,include_windspeed_output,include_stressx_output,include_stressy_output,include_Etide_output,include_Ewind_output,include_u_mean_surface_output,include_u_mean_bottom_output,include_grow1_mean_surface_output,include_grow1_mean_bottom_output,include_uptake1_mean_surface_output,include_uptake1_mean_bottom_output,include_tpn1_output,include_tpg1_output,include_speed3_output)
        # results,errors = pool.map(func, range(num_lines))
        results, errors = zip(*pool.map(func, range(len(lat_domain))))
        # results = pool.map(func, range(num_lines))

        if generate_netcdf_files:

            # run_start_date = str(year)+'-01-01'
            # df = pd.DataFrame(columns=(column_names))
            # i=0
            # for result in results:
            #     lines = result.split('\n')[:-1]
            #     for line in lines:
            #         # print line
            #         df.loc[i] = map(float,line.split())
            #         i+=1
            #
            # for column_name in column_names[4::]:
            #     output_cube = put_data_into_cube(df,df_domain,column_name,specifying_names,standard_name,long_name,var_name,units,run_start_date)
            #     iris.fileformats.netcdf.save(output_cube, output_directory+output_file_name+'_'+column_name.replace(" ", "")+'_'+str(year)+'.nc', zlib=True, complevel=2)

            run_start_date = str(year)+'-01-01'
            df = pd.DataFrame(columns=(column_names))
            i=0
            tmp_array = np.zeros([len(column_names),np.sum([len(result.split('\n')[:-1]) for result in results])])
            for result in results:
                lines = result.split('\n')[:-1]
                for line in lines:
                    # print line
                    # df.loc[i] = map(float,line.split())
                    tmp_array[:,i] = map(float,line.split())
                    i+=1

#             df = pd.DataFrame({column_names[0]: tmp_array[0,:], column_names[1]: tmp_array[1,:], column_names[2]: tmp_array[2,:], column_names[3]: tmp_array[3,:], column_names[4]: tmp_array[4,:], column_names[5]: tmp_array[5,:], column_names[6]: tmp_array[6,:], column_names[7]: tmp_array[7,:], column_names[8]: tmp_array[8,:]})
# need to make this generic based on no column_names
            df = pd.DataFrame({column_names[0]: tmp_array[0,:]})

            for i in range(len(column_names)-1):
                df[column_names[i+1]] = tmp_array[i+1,:]

            # uncomment this when this set of runs is complete - needed for runs wich span the 0 lon line
            #df.longitude.values[np.where(df.longitude.values >= 180)] -= 360

            func = partial(output_netcdf,year,column_names,df,df_domain,specifying_names,standard_name,long_name,var_name,units,run_start_date, output_directory,output_file_name)
            my_log = zip(*pool.map(func, range(4,len(column_names))))
        else:
            with open(output_directory+output_file_name+'_'+str(year),'w') as fout:
                for result in results:
                    fout.write(result)
            if write_error_output:
                with open(output_directory+output_file_name+'_error_'+str(year),'w') as fout:
                    for error in errors:
                        fout.write(error)
        pool.close()

    if not parallel_processing:
        # non parallel version
        with open(output_directory+output_file_name+'_'+str(year),'w') as fout:
            for i in range(len(lat_domain)):
                out,err = run_model(domain_file_name, lats_lons, year, start_year, unique_job_id, met_data_temporary_location,lon_domain,lat_domain,smaj1,smin1,smaj2,smin2,smaj3,smin3,smaj4,smin4,smaj5,smin5,woa_nutrient,alldepth,include_depth_output,include_temp_surface_output,include_temp_bottom_output,include_chlorophyll_surface_output,include_phyto_biomass_surface_output,include_phyto_biomass_bottom_output,include_PAR_surface_output,include_PAR_bottom_output,include_windspeed_output,include_stressx_output,include_stressy_output,include_Etide_output,include_Ewind_output,include_u_mean_surface_output,include_u_mean_bottom_output,include_grow1_mean_surface_output,include_grow1_mean_bottom_output,include_uptake1_mean_surface_output,include_uptake1_mean_bottom_output,include_tpn1_output,include_tpg1_output,include_speed3_output,i)
                fout.write(out)



    #clean up and leftover met files
    try:
        files_to_delete = glob.glob(met_data_temporary_location+'*.dat')
        [os.remove(f) for f in files_to_delete]
    except:
        print 'no met files to clean up'



remove_files = glob.glob(base_directory+'main/*'+unique_job_id+'*')
try:
    remove_files.remove(base_directory+'main/restart'+unique_job_id+'.dat')
except:
    pass
for remove_file in remove_files:
    os.remove(remove_file)
