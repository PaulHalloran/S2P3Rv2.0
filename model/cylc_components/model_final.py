import glob
import toml
import os
import pickle

def main():
    
    # Load config
    cylc_components_dir = "/net/home/h06/jwiggs/shelf_seas_regional/S2P3Rv2.0/model/cylc_components"
    config = toml.load(os.path.join(cylc_components_dir, "config.toml"))

    unique_job_id = pickle.load(open(os.path.join(cylc_components_dir, "run_id.p"), "rb"))

    remove_files = glob.glob(config["base_directory"] + "main/*" + unique_job_id + "*")
    try:
        remove_files.remove(config["base_directory"] + "main/restart" + unique_job_id + ".dat")
    except:
        pass
    for remove_file in remove_files:
        os.remove(remove_file)
