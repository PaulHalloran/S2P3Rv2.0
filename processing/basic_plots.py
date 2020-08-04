import iris
import glob
import matplotlib.pyplot as plt
import iris.plot as iplt
import iris.coord_categorisation
import numpy as np

#####
# Edit the three lines bellow to specify model output files. The example here is for files like: era5_uk_surfacetemperature_2000.nc, era5_uk_surfacetemperature_2001.nc, era5_uk_surfacetemperature_2002.nc
###

directory = './'
file_name = 'era5_uk'
variable = 'surfacetemperature'

#####
# Functions used to read in and plot the model data
####

def read_in_model_data(directory,file_name,variable):
    # Note that reading the data in with Iris is not completely trivial because of the time units of the files
    files = glob.glob(directory+'/'+file_name+'_'+variable+'_????.nc')
    files.sort()
    cubes = iris.load(files)
    #unifying time coordinate
    [cube.coord('time').convert_units('days since 1800-01-01 00:00:0.0') for cube in cubes]
    cube = cubes.concatenate_cube()
    return cube

def plot_map_timemean(cube,mymin=False,mymax=False,output_filename=False):
    # Script to plot a contour map of the values averaged over the full time dimension of the supplied cube
    # Specify mymin and mymax to set a color range
    # Specify output_filename as the directory path and filename to save the plot. End in '.png' to save as a png file, or '.pdf' to save as a pdf
    if not((mymin == False) and (mymax == False)):
        iplt.pcolormesh(cube.collapsed('time',iris.analysis.MEAN),vmin=mymin,vmax=mymax)
    else:
        iplt.pcolormesh(cube.collapsed('time',iris.analysis.MEAN))
    if not(output_filename == False):
        plt.savefig(output_filename)
    plt.gca().coastlines('50m')
    plt.colorbar()
    plt.show()

def extract_range_of_years(cube,min_year,max_year):
    # Extracts data from within a range of years, inclusive of those specified years
    try:
        iris.coord_categorisation.add_year(cube,'time', name='year')
        # note that similar commandds are available for month, season etc.
    except:
        pass
    years = cube.coord('year').points
    mask = np.where((years >= min_year) & (years <= max_year))
    return cube[mask]


def extract_geographical_area(cube,lon_west,lon_east,lat_south,lat_north):
    # Extracts data within specified lat/lon bounds
    region = iris.Constraint(longitude=lambda v: lon_west <= v <=lon_east,latitude=lambda v: lat_south <= v <= lat_north)
    return cube.extract(region)



def area_average_to_timeseries(cube):
    # Apply an area weighted average across all data in the supplied region (over lat and lon)
    try:
        cube.coord('latitude').guess_bounds()
    except:
        pass
    try:
        cube.coord('longitude').guess_bounds()
    except:
        pass
    grid_areas = iris.analysis.cartography.area_weights(cube)
    return cube.collapsed(['longitude', 'latitude'],iris.analysis.MEAN, weights=grid_areas)


def plot_a_timeseries(cube,output_filename=False):
    # plot an area averaged cube as a timeseries
    # Specify output_filename as the directory path and filename to save the plot. End in '.png' to save as a png file, or '.pdf' to save as a pdf
    iplt.plot(cube)
    plt.xlabel('date')
    plt.ylabel('variable name')
    if not(output_filename == False):
        plt.savefig(output_filename)
    plt.show()


def average_daily_data_to_annual(cube):
    try:
        iris.coord_categorisation.add_year(cube, 'time', name='year')
    except:
        pass
    return cube.aggregated_by('year', iris.analysis.MEAN)

def average_daily_data_to_seasonal(cube):
    try:
        iris.coord_categorisation.add_season_year(cube, 'time', name='season_year')
    except:
        pass
    return cube.aggregated_by(['season_year','year'], iris.analysis.MEAN)

def average_daily_data_to_monthly(cube):
    try:
        iris.coord_categorisation.add_month(cube, 'time', name='month')
    except:
        pass
    return cube.aggregated_by(['month','year'], iris.analysis.MEAN)



#####
# Reading in the data
#####

cube = read_in_model_data(directory,file_name,variable)

#####
# Plotting the data (examples)
#####

# Plot a contourmap of all data averaged along the time dimension
plot_map_timemean(cube)

# Plot a contourmap of the data falling between 2000 and 2005 inclusive averaged along the time dimension
cube_2000_2005 = extract_range_of_years(cube,2000,2005)
plot_map_timemean(cube_2000_2005)

# Plot a contourmap of the data from 2000-2005 falling between the specified bounds, averaged along the time dimension
lon_west = 0.0
lon_east = 10.5
lat_south = 55.0
lat_north = 65.0
cube_2000_2005_region = extract_geographical_area(cube_2000_2005,lon_west,lon_east,lat_south,lat_north)
plot_map_timemean(cube_2000_2005_region)

# Extract a region then create a timeseries from mapped data by performing a weighted area average, then plot
lon_west = 0.0
lon_east = 10.5
lat_south = 55.0
lat_north = 65.0
cube_region = extract_geographical_area(cube,lon_west,lon_east,lat_south,lat_north)
cube_timeseries = area_average_to_timeseries(cube_region)
plot_a_timeseries(cube_timeseries)


# plot daily, monthly, seasonal and annual mean cube_timeseries
# Annual
cube_timeseries = area_average_to_timeseries(cube_region)
cube_timeseries_annual = average_daily_data_to_annual(cube_timeseries)
plot_a_timeseries(cube_timeseries_annual)

# Seasonal
cube_timeseries_seasonal = average_daily_data_to_seasonal(cube_timeseries)
plot_a_timeseries(cube_timeseries_seasonal)

# Monthly
cube_timeseries_monthly = average_daily_data_to_monthly(cube_timeseries)
plot_a_timeseries(cube_timeseries_monthly)
