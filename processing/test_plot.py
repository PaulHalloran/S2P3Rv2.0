import iris
import glob
import matplotlib.pyplot as plt
import iris.plot as iplt
import iris.coord_categorisation
import numpy as np

def read_in_model_data(directory,file_name,variable):
    # Note that reading the data in with Iris is not completely trivial because of the time units of the files
    files = glob.glob(directory+'/'+file_name+'_'+variable+'_????.nc')
    files.sort()
    cubes = iris.load(files)
    #unifying time coordinate
    [cube.coord('time').convert_units('days since 1800-01-01 00:00:0.0') for cube in cubes]
    if len(cubes) > 1:
        cube = cubes.concatenate_cube()
    else:
        cube = cubes[0]
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


directory = '../test_output/'
file_name = 'test_output'
variable = 'surfacetemperature'

cube = read_in_model_data(directory,file_name,variable)

plot_map_timemean(cube)

