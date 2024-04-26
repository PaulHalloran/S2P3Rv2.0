import argparse
import glob
import multiprocessing as mp
import os
import shutil
import subprocess
import sys
import time
import uuid
import pickle
from functools import partial
from math import asin, cos, sqrt

import pandas as pd

import toml
from model_utils import output_netcdf, put_data_into_cube, run_model

def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "--year",
        type=int,
        help="Year to run",
    )
    args = parser.parse_args(args)
    return args

##################################################
# you may need to change things here             #
##################################################

cylc_components_dir = "/net/home/h06/jwiggs/shelf_seas_regional/S2P3Rv2.0/model/cylc_components"
config = toml.load(os.path.join(cylc_components_dir, "config.toml"))

args = parse_args(sys.argv[1:])

num_procs = (
    mp.cpu_count()
)  # this will use all available processors. Note that on a multi-node machine it can only use the processors on one node
# num_procs = 1 # The default is to use all available processors, but it is possible to specify the number of processors.

columns = [
    config['include_depth_output'],
    config['include_temp_surface_output'],
    config['include_temp_bottom_output'],
    config['include_chlorophyll_surface_output'],
    config['include_phyto_biomass_surface_output'],
    config['include_phyto_biomass_bottom_output'],
    config['include_PAR_surface_output'],
    config['include_PAR_bottom_output'],
    config['include_windspeed_output'],
    config['include_stressx_output'],
    config['include_stressy_output'],
    config['include_Etide_output'],
    config['include_Ewind_output'],
    config['include_u_mean_surface_output'],
    config['include_u_mean_bottom_output'],
    config['include_grow1_mean_surface_output'],
    config['include_grow1_mean_bottom_output'],
    config['include_uptake1_mean_surface_output'],
    config['include_uptake1_mean_bottom_output'],
    config['include_tpn1_output'],
    config['include_tpg1_output'],
    config['include_speed3_output'],
]
column_names_all = [
    "depth",
    "surface temperature",
    "bottom temperature",
    "surface chlorophyll",
    "surface phyto biomass",
    "bottom phyto biomass",
    "surface PAR",
    "bottom PAR",
    "windspeed",
    "stress_x",
    "stress_y",
    "Etide",
    "Ewind",
    "u_mean_surface",
    "u_mean_bottom",
    "grow1_mean_surface",
    "grow1_mean_bottom",
    "uptake1_mean_surface",
    "uptake1_mean_bottom",
    "tpn1",
    "tpg1",
    "speed3",
]

if config['generate_netcdf_files']:
    import numpy as np
    import iris
    import pandas as pd
    from itertools import compress
    from cf_units import Unit
    from iris import coords as iris_coords
    from iris import cube as iris_cube

    column_names = ["day", "longitude", "latitude"] + list(
        compress(column_names_all, map(bool, columns))
    )
    specifying_names = False
    ## If specifying_names above is set to True, specify the below. If not, ignore ##
    standard_name = [
        "sea_surface_temperature",
        "sea_surface_temperature",
        "sea_surface_temperature",
        "sea_surface_temperature",
        "sea_surface_temperature",
    ]
    long_name = [
        "Sea Surface Temperature",
        "Sea Surface Temperature",
        "Sea Surface Temperature",
        "Sea Surface Temperature",
        "Sea Surface Temperature",
    ]
    var_name = ["tos", "tos", "tos", "tos", "tos"]
    units = ["K", "K", "K", "K", "K"]
    if not (specifying_names):
        standard_name = np.tile("sea_surface_temperature", len(column_names))
        long_name = np.tile("Sea Surface Temperature", len(column_names))
        var_name = np.tile("tos", len(column_names))
        units = np.tile("K", len(column_names))


# Load pickled variables
unique_job_id = pickle.load(open(os.path.join(cylc_components_dir, "run_id.p"), "rb"))
lats_lons = pickle.load(open(os.path.join(cylc_components_dir, "lats_lons.p"), "rb"))
df_domain = pickle.load(open(os.path.join(cylc_components_dir, "domain_df.p"), "rb"))
(
        lon_domain,
        lat_domain,
        alldepth,
        smaj1,
        smin1,
        smaj2,
        smin2,
        smaj3,
        smin3,
        smaj4,
        smin4,
        smaj5,
        smin5,
        woa_nutrient,
    ) = pickle.load(open(os.path.join(cylc_components_dir, "tide_nutrient_lists.p"), "rb"))

##################################################
# main program                                   #
##################################################

#num_lines = sum(1 for line in open(base_directory + "domain/" + domain_file_name)) - 1
# num_lines = 10

os.chdir(os.path.join(config["base_directory"], "main"))

# year = start_year
print(args.year)
# clean up and prexisting met files
try:
    files_to_delete = glob.glob(config["met_data_temporary_location"] + "*.dat")
    [os.remove(f) for f in files_to_delete]
except:
    print("no met files to clean up")

subprocess.call(
    "tar -C "
    + config["met_data_temporary_location"]
    + " -zxf "
    + config["met_data_location"]
    + "met_data_"
    + str(args.year)
    + ".tar.gz",
    shell=True,
)
try:
    shutil.move(
        config["output_directory"] + config["output_file_name"] + "_" + str(args.year),
        config["output_directory"] + config["output_file_name"] + "_" + str(args.year) + "_previous",
    )
except:
    print("no previous output text file to move")

for column_name in column_names:
    try:
        shutil.move(
            config["output_directory"]
            + config["output_file_name"]
            + "_"
            + column_name.replace(" ", "")
            + "_"
            + str(args.year)
            + ".nc",
            config["output_directory"]
            + config["output_file_name"]
            + "_"
            + column_name.replace(" ", "")
            + "_"
            + str(args.year)
            + "_previous"
            + ".nc",
        )
    except:
        print("no previous " + column_name + " output netcdf file to move")

if config["parallel_processing"]:
    pool = mp.Pool(processes=num_procs)
    func = partial(
        run_model,
        config['executable_file_name'],
        config["domain_file_name"],
        config['nutrient_file_name'],
        lats_lons,
        args.year,
        config["start_year"],
        unique_job_id,
        config["met_data_temporary_location"],
        lon_domain,
        lat_domain,
        smaj1,
        smin1,
        smaj2,
        smin2,
        smaj3,
        smin3,
        smaj4,
        smin4,
        smaj5,
        smin5,
        woa_nutrient,
        alldepth,
        config["include_depth_output"],
        config["include_temp_surface_output"],
        config["include_temp_bottom_output"],
        config["include_chlorophyll_surface_output"],
        config["include_phyto_biomass_surface_output"],
        config["include_phyto_biomass_bottom_output"],
        config["include_PAR_surface_output"],
        config["include_PAR_bottom_output"],
        config["include_windspeed_output"],
        config["include_stressx_output"],
        config["include_stressy_output"],
        config["include_Etide_output"],
        config["include_Ewind_output"],
        config["include_u_mean_surface_output"],
        config["include_u_mean_bottom_output"],
        config["include_grow1_mean_surface_output"],
        config["include_grow1_mean_bottom_output"],
        config["include_uptake1_mean_surface_output"],
        config["include_uptake1_mean_bottom_output"],
        config["include_tpn1_output"],
        config["include_tpg1_output"],
        config["include_speed3_output"],
    )
    # results,errors = pool.map(func, range(num_lines))
    results, errors = zip(*pool.map(func, range(len(lat_domain))))
    # results = pool.map(func, range(num_lines))

    # with open(output_directory+output_file_name+'_error_'+str(year),'w') as fout:
    #            for error in errors:
    #                fout.write(str(error))

    if config["generate_netcdf_files"]:

        # run_start_date = str(year)+'-01-01'
        # df = pd.DataFrame(columns=(column_names))
        # i=0
        # for result in results:
        #     lines = result.split('\n')[:-1]
        #     for line in lines:
        #         # print(line
        #         df.loc[i] = map(float,line.split())
        #         i+=1
        #
        # for column_name in column_names[4::]:
        #     output_cube = put_data_into_cube(df,df_domain,column_name,specifying_names,standard_name,long_name,var_name,units,run_start_date)
        #     iris.fileformats.netcdf.save(output_cube, output_directory+output_file_name+'_'+column_name.replace(" ", "")+'_'+str(year)+'.nc', zlib=True, complevel=2)

        run_start_date = str(args.year) + "-01-01"
        df = pd.DataFrame(columns=(column_names))
        i = 0
        tmp_array = np.zeros(
            [
                len(column_names),
                np.sum([len(result.split(b"\n")[:-1]) for result in results]),
            ]
        )
        for result in results:
            result = result.decode("ascii")
            lines = result.split("\n")[:-1]
            for line in lines:
                # print(line
                # df.loc[i] = map(float,line.split())
                tmp_array[:, i] = list(map(float, line.split()))
                i += 1

        #             df = pd.DataFrame({column_names[0]: tmp_array[0,:], column_names[1]: tmp_array[1,:], column_names[2]: tmp_array[2,:], column_names[3]: tmp_array[3,:], column_names[4]: tmp_array[4,:], column_names[5]: tmp_array[5,:], column_names[6]: tmp_array[6,:], column_names[7]: tmp_array[7,:], column_names[8]: tmp_array[8,:]})
        # need to make this generic based on no column_names
        df = pd.DataFrame({column_names[0]: tmp_array[0, :]})

        for i in range(len(column_names) - 1):
            df[column_names[i + 1]] = tmp_array[i + 1, :]

        # uncomment this when this set of runs is complete - needed for runs wich span the 0 lon line
        df.longitude.values[np.where(df.longitude.values >= 180)] -= 360

        func = partial(
            output_netcdf,
            args.year,
            column_names,
            df,
            df_domain,
            specifying_names,
            standard_name,
            long_name,
            var_name,
            units,
            run_start_date,
            config["output_directory"],
            config["output_file_name"],
        )
        my_log = zip(*pool.map(func, range(4, len(column_names))))
    else:
        with open(
            config["output_directory"] + config["output_file_name"] + "_" + str(args.year), "w"
        ) as fout:
            for result in results:
                fout.write(result)
        if config["write_error_output"]:
            with open(
                config["output_directory"] + config["output_file_name"] + "_error_" + str(args.year), "w"
            ) as fout:
                for error in errors:
                    fout.write(error)
    pool.close()

# clean up and leftover met files
try:
    files_to_delete = glob.glob(config["met_data_temporary_location"] + "*.dat")
    [os.remove(f) for f in files_to_delete]
except:
    print("no met files to clean up")


