#################
# edit the following three lines
#################

output_file_name = 'initial_nitrate.dat'
domain_file_name = 's12_m2_s2_n2_h_map.dat'
location_of_World_Ocean_Atlas13_Nitrate_file =  './woa13_all_n13_01.nc'
#woa13_all_n13_01.nc is the N. hemisphere winter file
#annual, summer or winter netcdf file (note, use winter file for better values in N. hem and summer file for better value in S. hem) from https://www.nodc.noaa.gov/cgi-bin/OC5/woa13/woa13oxnu.pl

import iris
import numpy as np
import pandas as pd
import tempfile
import shutil
import os
import subprocess
import csv
from math import cos, asin, sqrt
from scipy.spatial import KDTree
import math

def cartesian(latitude, longitude, elevation = 0):
    # Convert to radians
    latitude = latitude * (math.pi / 180)
    longitude = longitude * (math.pi / 180)
    R = 6371 # 6378137.0 + elevation  # relative to centre of the earth
    X = R * math.cos(latitude) * math.cos(longitude)
    Y = R * math.cos(latitude) * math.sin(longitude)
    Z = R * math.sin(latitude)
    return (X, Y, Z)

def find_loc(lat, lon):
    # cartesian_coord = cartesian(lat, lon)
    # closest = tree.query([cartesian_coord])
    # cartesian_coord = cartesian(lat, lon)
    closest = tree.query([lat, lon])
    return closest[1]


###########################
# read in and process the world ocean atlas nitrate
###########################

cube=iris.load(location_of_World_Ocean_Atlas13_Nitrate_file)[3][0][0:24]
#just taking the top 25 levels so we are only considering the top 200m of the water column

cube2 = cube.copy()
cube2 = cube2[0]
cube2.data.data[:]=cube2.data.fill_value

new_data = cube2.data

print 'constructing bottom water nitrate array'
for i in range(np.shape(cube2)[0]):
    print 'processing latitude ',i
    for j in range(np.shape(cube2)[1]):
        try:
            new_data[i,j] = cube.data[cube.data.mask[:,i,j]==False][-1,i,j]
        except:
            pass


cube2.data = new_data

#qplt.contourf(cube2,31)
#plt.savefig('woa_bottom_N.png')


##################################
# extract data for the chosen domain
##################################

# tides_df = pd.read_csv(domain_file_name,names=['lon','lat','t1','t2','t3','t4','t5','t6','t7','t8','t9','t10','depth'],delim_whitespace=True,skiprows=1)
fwidths=[8,8,6,6,6,6,6,6,6,6,6,6,8]
print 'reading in lats and lons from domain file'
tides_df = pd.read_fwf(domain_file_name,names=['lon','lat','t1','t2','t3','t4','t5','t6','t7','t8','t9','t10','depth'],widths = fwidths,
                 skiprows=[0],dtype={'lon':float,'lat':float,'t1':float,'t2':float,'t3':float,'t4':float,'t5':float,'t6':float,'t7':float,'t8':float,'t9':float,'t10':float,'depth':float},usecols=['lon','lat','depth'])


#if min(tides_df.lon) < 0.0:
#    tides_df.lon = 360.0+tides_df.lon

lon_west = np.min(tides_df['lon'].values)
lon_east = np.max(tides_df['lon'].values)
lat_south = np.min(tides_df['lat'].values)
lat_north = np.max(tides_df['lat'].values)

cube_region_tmp = cube2.intersection(longitude=(lon_west, lon_east))
cube2 = cube_region_tmp.intersection(latitude=(lat_south, lat_north))

nitrate = np.zeros(len(tides_df['lon'].values))
nitrate[:] = np.nan

cube_lats = cube2.coord('latitude').points
cube_lons = cube2.coord('longitude').points
cube_data = cube2.data

cube_lats_b,cube_lons_b  = np.meshgrid(cube_lats, cube_lons, indexing='ij')

cube_data2 = np.reshape(cube_data,np.shape(cube_data)[0]*np.shape(cube_data)[1])
cube_lats2 = np.reshape(cube_lats_b,np.shape(cube_data)[0]*np.shape(cube_data)[1])
cube_lons2 = np.reshape(cube_lons_b,np.shape(cube_data)[0]*np.shape(cube_data)[1])

cube2.data = np.ma.masked_array(cube2.data)

places = []
for i in range(len(cube2.coord('longitude').points)):
    for j in range(len(cube2.coord('latitude').points)):
        if not cube2.data.mask[j,i]:
            coordinates = [cube2.coord('latitude').points[j], cube2.coord('longitude').points[i]]
            # cartesian_coord = cartesian(*coordinates)
            # places.append(cartesian_coord)
            places.append(coordinates)


tree = KDTree(places)

for i in range(len(tides_df['lon'].values)):
#    print i,' of ',len(tides_df['lon'].values)
    lon = tides_df['lon'].values[i]
    lat = tides_df['lat'].values[i]
    location = places[find_loc(lat, lon)]
    index_lon = cube2.coord('longitude').nearest_neighbour_index(location[1])
    index_lat = cube2.coord('latitude').nearest_neighbour_index(location[0])
    nitrate[i] = cube2.data[index_lat,index_lon]


##################################
# write the new initialisation file
##################################


output_df = tides_df.copy()
output_df['nitrate'] = nitrate

try:
    os.remove(output_file_name)
except:
    print ('no file to remove')

use_cols=['lon','lat','nitrate']

outF = open(output_file_name, "w")
# outF.write('1\n')
for j in range(len(output_df)):
   format=["%8.3f","%8.3f","%6.1f"]
   outF.write(''.join(format[i] % output_df[use_cols[i]].iloc[j] for i in range(3)))
  # write line to output file
   outF.write("\n")
outF.close()
