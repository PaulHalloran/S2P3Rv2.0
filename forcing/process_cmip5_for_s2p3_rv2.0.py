
##################################
# importing modules we will be using later
##################################

import iris
import numpy as np
import matplotlib.pyplot as plt
import glob
import pandas as pd
import subprocess
import os
import csv
from math import cos, asin, sqrt
import threading
import time
import iris.coord_categorisation
from scipy import interpolate
from time import sleep
import sys
from scipy.interpolate import RectBivariateSpline
import scipy.interpolate as interp
import multiprocessing as mp
from functools import partial
import tarfile
from scipy.spatial import KDTree

##################################
# things that may need to be edited
##################################


min_depth_lim = 10.0
max_depth_lim = 100.0

start_year = 2006
end_year = 2007

cmip_model = 'HadGEM2-ES' # note that this should match exactly the model name used in the filename
experiment = 'rcp85' # note that this shoudl match exactly the experiment name used in the filename

# SOME LON LATS DON'T HAVE METEROLOGY: /data/NAS-ph290/ph290/cmip5/for_s2p3_rv2.0/merged/
# This appears to be at the edge of the dataset... maybe extract a larger regoin than is required


domain_file = 's12_m2_s2_n2_h_map.dat'
# Note, the script will fail with the error 'OverflowError: cannot serialize a string larger than 2 GiB'
# if this file is too big. For example, a global 4km dataste is too big, but 4km from 30S to 30N is OK
# The limitation is that the data has to be pickled to be run in a parallellised way, and currently
# pickle has a 2GB limit in python2

#Specify where the ncep data is stored on your computer
#directory_containing_files_to_process = '/data/NAS-ph290/ph290/cmip5/for_s2p3_rv2.0/merged/'
directory_containing_files_to_process = '/data/BatCaveNAS/ph290/tmp/'
output_directory = '/home/ph290/Downloads/'

#value to set the minimum wind value to to avoid build up of heat in high cloud, low wind speed situations
# min_wind_value = 2.0
# just for ref. model requires: idum,wind_speed,wind_dir,cloud,airT,airP,humid

##################################
# some functions we will make use of later
##################################

def progress_bar(k,n):
    #k=current itteration, n= number of itterations
    sys.stdout.write("\r%d%%" % int((float(k)/n)*100.0))
    sys.stdout.flush()


def land_fill(mask_cube,cubes):
    ########## set land point mask values to True ##########
    #I think the below will make sure that land is masked for variables provided on the a c-grid, e.g. u and v in HadGEM2
    # http://cdn.intechopen.com/pdfs/43438/InTech-Grids_in_numerical_weather_and_climate_models.pdf
    #PRESENTLY UNTESTED!!!
    if np.shape(cubes)[0] == np.shape(mask_cube)[0]-1:
        data = mask_cube.data.copy()
        data = data + np.roll(data,1,axis=0)
        mask_cube = cubes[0].copy()
        mask_cube.data = data[0:-1][:]
    where_mask = mask_cube.data > 0.0
    where_mask_3d = np.repeat(where_mask[np.newaxis, :, :], np.shape(cubes[0])[0], axis=0)
    data = []
    for i,cube in enumerate(cubes):
        cubes[i].data.mask = where_mask_3d
        data.append(cubes[i].data)
    ########## filling land points with nearest neighbour value ##########
    for i,dummy1 in enumerate(data):
        print 'processing ',i,' out of ',len(data),' variables'
        a = data[i][0,:,:].copy()
        x,y=np.mgrid[0:a.shape[0],0:a.shape[1]]
        xygood = np.array((x[~a.mask],y[~a.mask])).T
        xybad = np.array((x[a.mask],y[a.mask])).T
        nearest = KDTree(xygood).query(xybad)[1]
        for j in range(data[i].shape[0]):
            # print 'filling ',j,' out of ',data[0].shape[0],' 2d fields'
            data[i][j,a.mask] = data[i][j,~a.mask][nearest]
        cubes[i].data = data[i]
        cubes[i].data.mask[:] = False
    return cubes


# def huss_to_hurs_old(huss,tas,psl):
#     #an approximation...
#     # https://github.com/pwcazenave/fvcom-toolbox/blob/master/utilities/qair2rh.m
#     temp = tas - 273.15
#     es = 6.112 * np.exp((17.67 * temp) / (temp + 243.5))
#     e = huss * (psl * 0.01) / (0.378 * huss + 0.622)
#     rh = e / es
#     rh[np.where(rh > 1.0)] = 1.0
#     rh[np.where(rh < 0.0)] = 0.0
#     return rh

def huss_to_hurs(huss,tas,psl):
    t = tas-273.15
    # t = 20.0
    # psl = 1011.0
    # huss=0.002
    sat_vap=10.0**((0.7859+0.03477*t)/(1.0+0.00412*t)) # this is consistant with teh definition in the model
    vap = (huss/0.622) * psl
    rh = vap/sat_vap
    return rh


##################################
# Other things that need to be defined, but probably not changed
##################################

output_filename = 'meterological_data'

if len(glob.glob(output_directory+output_filename+'*.dat')) <> 0:
    print 'files already exist in the output directory. Please move, delete or point output to a new directory'
    print 'first file'
    print glob.glob(output_directory+output_filename+'*.dat')[0]
    sys.exit()

if len(glob.glob(output_directory+'met_data_*.tar.gz')) <> 0:
    print 'files already exist in the output directory. Please move, delete or point output to a new directory'
    print 'first file'
    print glob.glob(output_directory+'met_data_*.tar.gz')[0]
    sys.exit()

# df = pd.read_csv(domain_file,names=['lon','lat','t1','t2','t3','t4','t5','t6','t7','t8','t9','t10','depth'],delim_whitespace=True,skiprows=[0],dtype={'lon':float,'lat':float,'t1':float,'t2':float,'t3':float,'t4':float,'t5':float,'t6':float,'t7':float,'t8':float,'t9':float,'t10':float,'depth':float})
fwidths=[8,8,6,6,6,6,6,6,6,6,6,6,8]
print 'reading in lats and lons from domain file'
df = pd.read_fwf(domain_file,names=['lon','lat','t1','t2','t3','t4','t5','t6','t7','t8','t9','t10','depth'],widths = fwidths,
                 skiprows=[0],dtype={'lon':float,'lat':float,'t1':float,'t2':float,'t3':float,'t4':float,'t5':float,'t6':float,'t7':float,'t8':float,'t9':float,'t10':float,'depth':float},usecols=['lon','lat','depth'])

print 'completed reading in lats and lons from domain file'

# input_variables = ['vwnd','uwnd','tcdc','rhum','air.2m','pres'] #these are the names of the variables in the NCEP files
input_variables = ['vas','uas','clt','huss','tas','psl','wind_speed'] #these are the cmip5 names, and should all be avaiable daily NOTE WE HAVE huss not hurs... so need to convert huss to hurs
# North/South wind vector, East/West wind vector, Total cloud cover, relative humidity, 2m air temperature, sea level pressure
#conversion specific to relative: http://www.whoi.edu/page.do?pid=30578

output_filename = 'meterological_data'

##################################
# pre-processing, reading, and storing sensibly the data from each of meterology variables, extracting data for the location of interest
##################################


sample_points_lat_lon = df.loc[(df['depth'] <= max_depth_lim) & (df['depth'] >= min_depth_lim)][['lat','lon']]
sample_points_lat_lon.drop_duplicates()

input_variables2 = np.array(input_variables).copy()
input_variables2 = np.append(input_variables2,'hurs')
# input_variables2 = np.append(input_variables2,'wind_speed')
input_variables2 = np.append(input_variables2,'wind_direction')

cube = iris.load_cube(directory_containing_files_to_process + cmip_model + '_' + input_variables[0]+'_'+experiment+'_all.nc')
iris.coord_categorisation.add_year(cube, 'time', name='year')
cube_year = cube[np.where(cube.coord('year').points == start_year)]
znew = np.zeros([np.shape(cube_year)[0],len(input_variables2),len(sample_points_lat_lon['lat'].values)])
znew_tmp = np.zeros([np.shape(cube_year)[0],len(sample_points_lat_lon['lat'].values)])

tmp = sample_points_lat_lon['lon'].values
tmp[np.where(tmp < 0.0)] = 360.0 + tmp[np.where(tmp < 0.0)]
sample_points_lat_lon['lon'] = tmp
# X,Y=np.meshgrid(cube_year.coord('latitude').points,cube_year.coord('longitude').points)



def interpolate_forcing_data(input_variables,sample_points_lat_lon,znew_tmp,cubes,k):
    single_input_variable = input_variables[k]
    cube_year = cubes[k]
    znew_tmp[:] = np.nan
    print 'processing '+single_input_variable
    for j in range(np.shape(cube_year)[0]):
        # met_cubes[single_input_variable]['day']={}
        # progress_bar(j,np.shape(cube_year)[0])
        # different interpolation options commented out
        # f = interpolate.interp2d(cube_year.coord('longitude').points, cube_year.coord('latitude').points, cube_year[j].data, kind='linear')
        # znew[j,k,:] = [f(sample_points_lat_lon['lon'].values[i], sample_points_lat_lon['lat'].values[i])[0] for i in range(len(sample_points_lat_lon['lat'].values))]
        #https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.RectBivariateSpline.html
        f = RectBivariateSpline(cube_year.coord('latitude').points,cube_year.coord('longitude').points, cube_year[j].data)
        znew_tmp[j,:] = [f(sample_points_lat_lon['lat'].values[i], sample_points_lat_lon['lon'].values[i])[0] for i in range(len(sample_points_lat_lon['lat'].values))]
        #https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.Rbf.html
        # f = interp.Rbf(X,Y, cube_year[j].data, function='linear', smooth=0)  # default smooth=0 for interpolation
        # znew[j,k,:] = [f(sample_points_lat_lon['lat'].values[i], sample_points_lat_lon['lon'].values[i]) for i in range(len(sample_points_lat_lon['lat'].values))]
    # return [znew_tmp.reshape(j*len(sample_points_lat_lon['lat'].values))]
    return [[znew_tmp]]


cwd = os.getcwd()

def ws_data_func(u_data, v_data):
    return np.sqrt( u_data**2 + v_data**2 )

def ws_units_func(u_cube, v_cube):
    if u_cube.units != getattr(v_cube, 'units', u_cube.units):
        raise ValueError("units do not match")
    return u_cube.units

# test if windspeed file_exists, if not produce it
exists = os.path.isfile(directory_containing_files_to_process + cmip_model + '_wind_speed_'+experiment+'_all.nc')
if not exists:
    u_cube = iris.load_cube(directory_containing_files_to_process + cmip_model + '_uas_'+experiment+'_all.nc')
    v_cube = iris.load_cube(directory_containing_files_to_process + cmip_model + '_vas_'+experiment+'_all.nc')
    ws_ifunc = iris.analysis.maths.IFunc(ws_data_func,ws_units_func)
    ws_cube = ws_ifunc(u_cube, v_cube, new_name='wind speed')
    iris.save(ws_cube, directory_containing_files_to_process + cmip_model + '_wind_speed_'+experiment+'_all.nc')

for year in range(start_year,end_year+1):
    print 'processing year ',year
    cubes = []
    cube_data=[]
    for k in range(len(input_variables)):
        single_input_variable = input_variables[k]
        print 'loading data for '+single_input_variable
        cube = iris.load_cube(directory_containing_files_to_process + cmip_model + '_' + single_input_variable+'_'+experiment+'_all.nc')
        iris.coord_categorisation.add_year(cube, 'time', name='year')
        cube_year = cube[np.where(cube.coord('year').points == year)]
        cube_data.append(cube_year.data)
        cubes.append(cube_year)
    #fill land grid boxes with values from nearest ocean grid box to avoid (e.g.) anomalously low winds speeds in some coastal grid boxes.
    mask_files = glob.glob(directory_containing_files_to_process + 'sftlf_fx_'+cmip_model+'_*.nc')
    if len(mask_files) == 0:
        print 'missing land-sea fraction file, variable sftlf, can not replace under land points with nearest under sea point'
    else:
        mask_file = mask_files[0] # select just one file if more have been downloaded
        mask_cube = iris.load_cube(mask_file)
        cubes = land_fill(mask_cube,cubes)
    znew[:] = np.NAN
    num_procs = mp.cpu_count()
    pool = mp.Pool(processes = np.min([6,num_procs]))
    func = partial(interpolate_forcing_data, input_variables,sample_points_lat_lon,znew_tmp,cubes)
    results = pool.map(func, range(len(input_variables)))
    znew[:,0:len(input_variables),:] = np.moveaxis(np.array(results)[:,0,0,:,:],1,0)
    # for k,single_input_variable in enumerate(input_variables):
        # cube = iris.load_cube(directory_containing_files_to_process + cmip_model + '_' + single_input_variable+'_'+experiment+'_all.nc')
        # iris.coord_categorisation.add_year(cube, 'time', name='year')
        # cube_year = cube[np.where(cube.coord('year').points == year)]
        # print 'processing days of year'
        # for j in range(np.shape(cube_year)[0]):
        #     # met_cubes[single_input_variable]['day']={}
        #     progress_bar(j,np.shape(cube_year)[0])
        #     # different interpolation options commented out
        #     # f = interpolate.interp2d(cube_year.coord('longitude').points, cube_year.coord('latitude').points, cube_year[j].data, kind='linear')
        #     #znew[j,k,:] = [f(sample_points_lat_lon['lon'].values[i], sample_points_lat_lon['lat'].values[i])[0] for i in range(len(sample_points_lat_lon['lat'].values))]
        #     #https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.RectBivariateSpline.html
        #     f = RectBivariateSpline(cube_year.coord('latitude').points,cube_year.coord('longitude').points, cube_year[j].data)
        #     znew[j,k,:] = [f(sample_points_lat_lon['lat'].values[i], sample_points_lat_lon['lon'].values[i])[0] for i in range(len(sample_points_lat_lon['lat'].values))]
        #     #https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.Rbf.html
        #     # f = interp.Rbf(X,Y, cube_year[j].data, function='linear', smooth=0)  # default smooth=0 for interpolation
        #     # znew[j,k,:] = [f(sample_points_lat_lon['lat'].values[i], sample_points_lat_lon['lon'].values[i]) for i in range(len(sample_points_lat_lon['lat'].values))]
    print 'calculating hurs'
    znew[:,np.where(input_variables2 == 'hurs')[0],:] = huss_to_hurs(znew[:,np.where(input_variables2 == 'huss')[0],:],znew[:,np.where(input_variables2 == 'tas')[0],:],znew[:,np.where(input_variables2 == 'psl')[0],:])
    #this simply converts the data in the 2m air temperature column from K to degrees C
    znew[:,np.where(input_variables2 == 'tas')[0],:] -= 273.15
    #this gets the pressure in the units used by the s2p3 model
    znew[:,np.where(input_variables2 == 'psl')[0],:] /= 100.0
    # relative humidity to percent from fraction
    # znew[:,np.where(input_variables2 == 'hurs')[0],:] *= 100.0
    #wind speed, using pythagoras (square root of the sum of the squares of the x and y vector give teg lenth of the 3rd side of the triangle)
    # znew[:,np.where(input_variables2 == 'wind_speed')[0],:] = np.sqrt(np.square(znew[:,np.where(input_variables2 == 'uas')[0],:]) + np.square(znew[:,np.where(input_variables2 == 'vas')[0],:]))
    #wind direction calculated using the function arctan2, then converted from radians to degrees
    znew[:,np.where(input_variables2 == 'wind_direction')[0],:] = np.rad2deg((np.arctan2(znew[:,np.where(input_variables2 == 'uas')[0],:],znew[:,np.where(input_variables2 == 'vas')[0],:])) + np.pi)
    #setting winds speeds below 2m/s to 2m/s to avoid issues with low wind speeds from predominantly land grid points
    # tmp = znew[:,np.where(input_variables2 == 'wind_speed')[0],:]
    # tmp2 = np.where(tmp < min_wind_value)
    # if not tmp2[0].size == 0:
    #     tmp[tmp2] = min_wind_value
    #     znew[:,np.where(input_variables2 == 'wind_speed')[0],:] = tmp
    print 'writing met data out'
    for u,longitude_point in enumerate(sample_points_lat_lon['lon'].values):
        latitude_point = sample_points_lat_lon['lat'].values[u]
        progress_bar(u,len(sample_points_lat_lon['lon'].values))
        # for v,latitude_point in enumerate(sample_points_lat_lon['lat'].values[0:500]):
        # delete.append(str(np.round(latitude_point,4))+str(np.round(longitude_point,4)))
        # delete2.append(str(int(latitude_point*10000))+str(int(longitude_point*10000)))
        # process the read-in data and get in the right format for output
        #To make things neater, we putting the data, which has so far been stored in a dictionary into a pandas dataframe
        #Once in this format it just makes other thinsg we may want to do easier. See https://www.tutorialspoint.com/python_pandas/python_pandas_dataframe.htm
        df2 = pd.DataFrame(data=znew[:,:,u],columns = input_variables2)
        #Essentially range(1, len(df) + 1) makes a list of numbers from 1 to x, where x is the length of the dataframe (the len(df) bit), i.e. the length of the meterological data we have
        #the ''[format(x, ' 5d') for x in r...' bit just takes each of those numbers one by one and formats them so istead of bineg '1', '2'... '10' etc. then are '    1', '    2'... '   10' etc. so that the colums line up correctly in the output file
        #This data is then stored in a new column in the dataframe, called 'day_number'
        df2['day_number'] = [format(x, ' 5d') for x in range(1, len(df2) + 1)]
        #The met file requires data is all to two decimal places, so the line below rounds all of the data to two decimal places.
        df2 = df2.round(2)
        ##################################
        # Write the data out to the file
        ##################################
        #this line simply writes out the olumns we are intersted in, in the order we are intersted in, in the firmat we are intersted in (2 decomal places, 10 characters between columns) to the file we specified at the start
        np.savetxt(output_directory+output_filename+'lat'+str(np.round(latitude_point,4))+'lon'+str(np.round(longitude_point,4))+'_'+str(year)+'.dat', df2[['day_number','wind_speed','wind_direction','clt','tas','psl','hurs']].values, fmt='%s%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f')
    pool.close()
    #tar and gzip the output files for each year:
    os.chdir(output_directory)
    names = [os.path.basename(x) for x in glob.glob(output_filename+'*.dat')]
    tar = tarfile.open('met_data_'+str(year)+'.tar.gz', 'w:gz')
    for name in names:
        tar.add(name)
    tar.close()
    #remove the files that have now been tar.gzped
    [os.remove(f) for f in names]
    os.chdir(cwd)
