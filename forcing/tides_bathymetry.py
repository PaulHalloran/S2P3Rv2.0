#note, latitudes go from -90 to 90
#minimum means furthest South
minimum_latitude = 10
maximum_latitude = 14
latitude_resolution = 0.1 #degrees

#note, longitudes go from -180 to 180
#minimum means furthest West
minimum_longitude = 90
maximum_longitude = 94
longitude_resolution = 0.1 #degrees

location_of_bathymetry_file =  './ETOPO1_Bed_g_gmt4.nc'

output_file_name = 's12_m2_s2_n2_h_map.dat'

import numpy as np
import pandas as pd
import tempfile
import shutil
import os
import subprocess
import csv
import iris


def replace_character_in_file(filename,character1,character2):
    with open(filename, 'r') as infile,open(filename+'cleaned_up', 'w') as outfile:
        data = infile.read()
        data = data.replace(character1, character2)
        outfile.write(data)
    subprocess.call(['mv '+filename+'cleaned_up'+' '+filename], shell=True)



def replace(file_path,file_path2, pattern, subst):
    #Create temp file
    fh, abs_path = tempfile.mkstemp()
    with os.fdopen(fh,'w') as new_file:
        with open(file_path) as old_file:
            for line in old_file:
                # new_file.write(line.replace(pattern[i], subst[i]))
                my_line = line
                for i,dummy in enumerate(pattern):
                    my_line = my_line.replace(pattern[i], subst[i])
                new_file.write(my_line)
    shutil.move(abs_path, file_path2)




def ap2ep(Au, PHIu, Av, PHIv):
    PHIu = PHIu / 180. * np.pi
    PHIv = PHIv / 180. * np.pi
    # Make complex amplitudes for u and v
    i = np.sqrt(-1+0j)
    u = Au * np.exp(-i * PHIu)
    v = Av * np.exp(-i * PHIv)
    # Calculate complex radius of anticlockwise and clockwise circles
    wp = (u + i * v) / 2. # for anticlockwise circles
    wm = np.conj(u - i * v) / 2. # for clockwise circles
    # and their amplitudes and angles
    Wp = np.abs(wp)
    Wm = np.abs(wm)
    THETAp = np.angle(wp)
    THETAm = np.angle(wm)
    # calculate the ellipse parameters
    SEMA = Wp + Wm
    SEMI = Wp - Wm
    ECC = SEMI / SEMA
    PHA = (THETAm - THETAp) / 2.
    INC = (THETAm + THETAp) / 2.
    PHA = PHA / np.pi * 180
    INC = INC / np.pi * 180
    THETAp = THETAp / np.pi * 180
    THETAm = THETAm / np.pi * 180
    THETAp[np.where(THETAp < 0)] = THETAp[np.where(THETAp < 0)] + 360
    THETAm[THETAm < 0] = THETAm[THETAm < 0] + 360
    PHA[PHA < 0] = PHA[PHA < 0] + 360
    INC[INC < 0] = INC[INC < 0] + 360
    return(PHA, INC, SEMI, SEMA)

tidal_components = ['m2', 's2', 'n2','k1','o1']
#options:m2,s2,n2,k2,k1,o1,p1,q1

# latitude_list = np.linspace(minimum_latitude,maximum_latitude,((maximum_latitude - minimum_latitude)/latitude_resolution)+latitude_resolution)
# longitude_list = np.linspace(minimum_longitude,maximum_longitude,(maximum_longitude - minimum_longitude)/longitude_resolution)
latitude_list = np.arange(minimum_latitude,maximum_latitude,latitude_resolution)
longitude_list = np.arange(minimum_longitude,maximum_longitude,longitude_resolution)

# if (minimum_latitude < 0.0) & (maximum_latitude > 0.0):
#     latitude_list = np.linspace(minimum_latitude,maximum_latitude,(maximum_latitude - minimum_latitude)/latitude_resolution+1.0)
#
# if (minimum_longitude < 0.0) & (maximum_longitude > 0.0):
#     longitude_list = np.linspace(minimum_longitude,maximum_longitude,(maximum_longitude - minimum_longitude)/longitude_resolution+1.0)


longitudes = []
latitudes = []

for lo in longitude_list:
    for la in latitude_list:
        longitudes.append(lo)
        latitudes.append(la)


data = {}
data['latitudes'] = np.array(latitudes).round(3)
data['longitudes'] = np.array(longitudes).round(3)
#Note that he time information below is completely arbitrary, and just included so that we can get the bathymetr data out of the predict executable
data['year'] = (np.zeros(np.size(latitudes))+2000).astype(int)
data['month'] = (np.zeros(np.size(latitudes)) + 10).astype(int)
data['day'] = (np.zeros(np.size(latitudes)) + 10).astype(int)
data['hour'] = (np.zeros(np.size(latitudes)) + 10).astype(int)
data['min'] = (np.zeros(np.size(latitudes)) + 10).astype(int)
data['second'] = (np.zeros(np.size(latitudes)) + 10).astype(int)

df = pd.DataFrame(data=data)

#this line simply writes out the columns we are intersted in, in teh order we are intersted in, in the firmat we are intersted in
# df[['latitudes','longitudes','year','month','day','hour','min','second']].to_csv('./lat_lon_time', index=False, header=False, float_format='%10.4f')
df[['latitudes','longitudes','year','month','day','hour','min','second']].to_csv('./lat_lon_time', index=False, header=False)
#unfortunately, I could not quickly find a way to successfully write the file without commas between the columns, so the line below simply strips the columns out from the columns.
replace_character_in_file('./lat_lon_time',',','   ')

output = {}

for tidal_component in tidal_components:
    replace('setup.inp_template','setup.inp', ['replace_this','and_swap_this','replace_tidal_constit'], ['u','output_1.out',tidal_component])
    subprocess.call(['./extract_HC<setup.inp'], shell=True)
    replace_character_in_file('output_1.out','       ************* Site is out of model grid OR land ***************','   0.000   0.000')
    replace('setup.inp_template','setup.inp', ['replace_this','and_swap_this','replace_tidal_constit'], ['v','output_2.out',tidal_component])
    subprocess.call(['./extract_HC<setup.inp'], shell=True)
    replace_character_in_file('output_2.out','       ************* Site is out of model grid OR land ***************','   0.000   0.000')
    # replace('setup.inp_template','setup.inp', ['replace_this','and_swap_this','replace_tidal_constit'], ['z','output_3.out',tidal_component])
    # subprocess.call(['./predict_tide<setup.inp'], shell=True)
    # replace_character_in_file('output_3.out','***** Site is out of model grid OR land *****','     10.10.2000 10:10:10     0.000     0.000')
    u_component = pd.read_csv('output_1.out', header=2, delimiter=r"\s+")
    v_component = pd.read_csv('output_2.out', header=2, delimiter=r"\s+")
    # bathy_component = pd.read_csv('output_3.out', header=3, delimiter=r"\s+")
    Au = u_component[tidal_component+'_amp'].values
    PHIu = u_component[tidal_component+'_ph'].values
    Av = v_component[tidal_component+'_amp'].values
    PHIv = v_component[tidal_component+'_ph'].values
    PHA, INC, SEMI, SEMA = ap2ep(Au, PHIu, Av, PHIv)
    output[tidal_component+'_SEMA'] = SEMA.round(1)
    output[tidal_component+'_SEMI'] = SEMI.round(1)


output['longitudes'] = df['longitudes'].round(3)
output['latitudes'] = df['latitudes'].round(3)
# output['depth'] = bathy_component['Depth(m)'].round(1)
# output['first_column'] = [" " for x in range(np.size(latitudes))]

cube = iris.load_cube(location_of_bathymetry_file)
cube_region_tmp = cube.intersection(longitude=(minimum_longitude, maximum_longitude))
cube2 = cube_region_tmp.intersection(latitude=(minimum_latitude, maximum_latitude))

latitude = iris.coords.DimCoord(np.arange(minimum_latitude, maximum_latitude+latitude_resolution, latitude_resolution), standard_name='latitude', units='degrees')
longitude = iris.coords.DimCoord(np.arange(minimum_longitude, maximum_longitude+longitude_resolution, longitude_resolution), standard_name='longitude', units='degrees')
regridding_cube = iris.cube.Cube(np.zeros((np.size(np.arange(minimum_latitude, maximum_latitude+latitude_resolution, latitude_resolution)), np.size(np.arange(minimum_longitude, maximum_longitude+longitude_resolution, longitude_resolution))), np.float32),standard_name='sea_surface_temperature', long_name='Sea Surface Temperature', var_name='tos', units='K',dim_coords_and_dims=[(latitude, 0), (longitude, 1)])

cube2_regridded = cube2.regrid(regridding_cube, iris.analysis.Nearest())

depths = np.zeros(len(output['latitudes']))
depths[:] = np.nan

# for i in range(np.shape(cube2_regridded)[1]):
#     for j in range(np.shape(cube2_regridded)[0]):
#         depths[(i * np.shape(cube2_regridded)[1]) + j] = cube2_regridded[j,i].data

lat = cube2_regridded.coord('latitude')
lon = cube2_regridded.coord('longitude')

for i in range(np.size(output['longitudes'])):
        lat_coord = lat.nearest_neighbour_index(output['latitudes'][i])
        lon_coord = lon.nearest_neighbour_index(output['longitudes'][i])
        depths[i] = cube2_regridded[lat_coord,lon_coord].data


depths = depths * -1.0

output['depth'] = depths

output_df = pd.DataFrame(data=output)
try:
    os.remove(output_file_name)
except:
    print ('no file to remove')

use_cols=['longitudes','latitudes','m2_SEMA','m2_SEMI','s2_SEMA','s2_SEMI','n2_SEMA','n2_SEMI','o1_SEMA','o1_SEMI','k1_SEMA','k1_SEMI','depth']

outF = open(output_file_name, "w")
outF.write('1\n')
#note removing -1 in for j in range(len(output_df)-1):
for j in range(len(output_df)):
   format=["%8.3f","%8.3f","%6.1f","%6.1f","%6.1f","%6.1f","%6.1f","%6.1f","%6.1f","%6.1f","%6.1f","%6.1f","%8.1f"]
   #note removing +1 in use_cols[i]].iloc[j+1]:
   outF.write(''.join(format[i] % output_df[use_cols[i]].iloc[j] for i in range(13)))
  # write line to output file
   outF.write("\n")
outF.close()
