import glob
import os
import pickle
import subprocess
import uuid

import pandas as pd

import toml


def main():

    # Load config
    cylc_components_dir = "/net/home/h06/jwiggs/shelf_seas_regional/S2P3Rv2.0/model/cylc_components"
    config = toml.load(os.path.join(cylc_components_dir, "config.toml"))

    # Remove old pickled objects
    if glob.glob(os.path.join(cylc_components_dir, "*.p")):
        for file in glob.glob(os.path.join(cylc_components_dir, "*.p")):
            os.remove(file)

    # Obtain filepath for domain and nurtient file
    domain_file_for_run = os.path.join(
        config["base_directory"], "domain", config["domain_file_name"]
    )
    nutrient_file_for_run = os.path.join(
        config["base_directory"], "domain", config["nutrient_file_name"]
    )

    # Create domain datafame
    df_domain = make_domain_df(domain_file_for_run)

    # Create tides and nutrients lists
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
    ) = make_tide_nutrient_lists(
        domain_file_for_run,
        nutrient_file_for_run,
        config["depth_min"],
        config["depth_max"],
    )

    # Make run id
    unique_job_id = generate_run_id()

    # Move met data for initial year
    extract_met_data(
        config["met_data_temporary_location"],
        config["met_data_location"],
        config["start_year"],
    )

    # Get list of initial met files
    initial_met_files = glob.glob(
        config["met_data_temporary_location"]
        + "*_"
        + str(config["start_year"])
        + ".dat"
    )
    lats_lons = make_lats_lons(initial_met_files)

    # Pickle required variables
    pickle_object(df_domain, os.path.join(cylc_components_dir, "domain_df.p"))
    pickle_object(
        [
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
        ],
        os.path.join(cylc_components_dir, "tide_nutrient_lists.p"),
    )
    pickle_object(unique_job_id, os.path.join(cylc_components_dir, "run_id.p"))
    pickle_object(lats_lons, os.path.join(cylc_components_dir, "lats_lons.p"))


def make_domain_df(domain_file_for_run):
    fwidths = [8, 8, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 8]
    df_domain = pd.read_fwf(
        domain_file_for_run,
        names=[
            "lon",
            "lat",
            "t1",
            "t2",
            "t3",
            "t4",
            "t5",
            "t6",
            "t7",
            "t8",
            "t9",
            "t10",
            "depth",
        ],
        widths=fwidths,
        skiprows=[0],
        dtype={
            "lon": float,
            "lat": float,
            "t1": float,
            "t2": float,
            "t3": float,
            "t4": float,
            "t5": float,
            "t6": float,
            "t7": float,
            "t8": float,
            "t9": float,
            "t10": float,
            "depth": float,
        },
        usecols=["lon", "lat", "depth"],
    )

    return df_domain


def make_tide_nutrient_lists(
    domain_file_for_run, nutrient_file_for_run, depth_min, depth_max
):
    f = open(domain_file_for_run)
    lines = f.readlines()
    f2 = open(nutrient_file_for_run)
    lines2 = f2.readlines()
    lat_domain = []
    lon_domain = []
    alldepth = []
    smaj1 = []
    smin1 = []
    smaj2 = []
    smin2 = []
    smaj3 = []
    smin3 = []
    smaj4 = []
    smin4 = []
    smaj5 = []
    smin5 = []
    woa_nutrient = []
    counter = 0
    for i, line in enumerate(lines[1::]):
        depth = float(line[77:84])
        if (depth >= depth_min) & (depth <= depth_max) & (depth > 0.0):
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

    return (
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
    )


def pickle_object(obj, fname):
    pickle.dump(obj, open(fname, "wb"))


def generate_run_id():
    return str(uuid.uuid4())


def extract_met_data(met_data_temporary_location, met_data_location, year):
    subprocess.call(
        "tar -C "
        + met_data_temporary_location
        + " -zxf "
        + met_data_location
        + "met_data_"
        + str(year)
        + ".tar.gz",
        shell=True,
    )


def make_lats_lons(files):
    w, h = 2, len(files)
    lats_lons = [[0 for x in range(w)] for y in range(h)]
    for i, file in enumerate(files):
        tmp = file.split("lat")[-1].split(".dat")[0].split("lon")
        lats_lons[i][0] = float(tmp[0])
        lats_lons[i][1] = float(tmp[1].split("_")[0])

    return lats_lons


if __name__ == "__main__":
    main()
