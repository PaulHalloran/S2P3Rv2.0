import numpy as np
import subprocess
from cf_units import Unit

import iris
from iris import coords as iris_coords
from iris import cube as iris_cube


def put_data_into_cube(
    df,
    df_domain,
    variable,
    specifying_names,
    standard_name,
    long_name,
    var_name,
    units,
    run_start_date,
):
    latitudes = np.unique(df_domain["lat"].values)
    longitudes = np.unique(df_domain["lon"].values)
    latitudes_run = np.unique(df["latitude"].values)
    longitudes_run = np.unique(df["longitude"].values)
    times = np.unique(df["day"].values)
    latitude = iris_coords.DimCoord(
        latitudes, standard_name="latitude", units="degrees"
    )
    longitude = iris_coords.DimCoord(
        longitudes, standard_name="longitude", units="degrees"
    )
    # time = iris.coords.DimCoord(times, standard_name='time', units='days')
    time = iris_coords.DimCoord(
        times,
        standard_name="time",
        units=Unit("days since " + run_start_date + " 00:00:0.0", calendar="gregorian"),
    )
    if specifying_names:
        cube = iris_cube.Cube(
            np.full((times.size, latitudes.size, longitudes.size), -999.99, np.float32),
            standard_name=standard_name,
            long_name=long_name,
            var_name=var_name,
            units=units,
            dim_coords_and_dims=[(time, 0), (latitude, 1), (longitude, 2)],
        )
    else:
        cube = iris_cube.Cube(
            np.full((times.size, latitudes.size, longitudes.size), -999.99, np.float32),
            standard_name=None,
            long_name=None,
            var_name=None,
            units=None,
            dim_coords_and_dims=[(time, 0), (latitude, 1), (longitude, 2)],
        )
    # Z,X,Y = np.meshgrid(cube.coord('time').points,cube.coord('longitude').points,cube.coord('latitude').points)
    data = cube.data.copy()
    # data[:] = -999.99
    days = np.unique(df["day"].values)
    # shape = [X.shape[0],X.shape[2]]
    for i, day in enumerate(days):
        df_tmp = df.loc[df["day"] == day]
        for j, lat in enumerate(df_tmp["latitude"].values):
            lon = df_tmp["longitude"].values[j]
            lat_loc = np.where(
                np.around(cube.coord("latitude").points, decimals=6)
                == np.around(lat, decimals=6)
            )[0][0]
            lon_loc = np.where(
                np.around(cube.coord("longitude").points, decimals=6)
                == np.around(lon, decimals=6)
            )[0][0]
            data[i, lat_loc, lon_loc] = df_tmp[variable].values[j]
    data_tmp = np.asarray(data.data)
    data = np.ma.masked_where((data_tmp < -999.9) & (data_tmp > -1000.0), data)
    # data = np.ma.masked_where(data.data == 0.0,data)
    cube.data = data
    cube.data.fill_value = -999.99
    cube.data.data[~(np.isfinite(cube.data.data))] = -999.99
    cube.data = np.ma.masked_where(cube.data == -999.99, cube.data)
    return cube

def output_netcdf(
    year,
    column_names,
    df,
    df_domain,
    specifying_names,
    standard_name,
    long_name,
    var_name,
    units,
    run_start_date,
    output_directory,
    output_file_name,
    i,
):
    from iris import fileformats as iris_fileformats

    column_name = column_names[i]
    output_cube = put_data_into_cube(
        df,
        df_domain,
        column_name,
        specifying_names,
        standard_name,
        long_name,
        var_name,
        units,
        run_start_date,
    )
    iris_fileformats.netcdf.save(
        output_cube,
        output_directory
        + output_file_name
        + "_"
        + column_name.replace(" ", "")
        + "_"
        + str(year)
        + ".nc",
        zlib=True,
        complevel=2,
    )
    return (
        output_directory
        + output_file_name
        + "_"
        + column_name.replace(" ", "")
        + "_"
        + str(year)
        + ".nc written"
    )

def run_model(
    executable_file_name,
    domain_file_name,
    nutrient_file_name,
    lats_lons,
    year,
    start_year,
    unique_job_id,
    met_data_temporary_location,
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
    include_depth_output,
    include_temp_surface_output,
    include_temp_bottom_output,
    include_chlorophyll_surface_output,
    include_phyto_biomass_surface_output,
    include_phyto_biomass_bottom_output,
    include_PAR_surface_output,
    include_PAR_bottom_output,
    include_windspeed_output,
    include_stressx_output,
    include_stressy_output,
    include_Etide_output,
    include_Ewind_output,
    include_u_mean_surface_output,
    include_u_mean_bottom_output,
    include_grow1_mean_surface_output,
    include_grow1_mean_bottom_output,
    include_uptake1_mean_surface_output,
    include_uptake1_mean_bottom_output,
    include_tpn1_output,
    include_tpg1_output,
    include_speed3_output,
    i,
):
    # modifying so that the fortran code looks for the correct met file, rather than us having to copy it into the working directory
    # lon,lat = return_domain_lon(base_directory+'domain/'+domain_file_name,i)
    lon_domain_tmp = float(lon_domain[i])
    if lon_domain_tmp < 0.0:
        lon_domain_tmp = 360.0 + lon_domain_tmp
    run_command = "\n".join(
        [
            "./{} << EOF".format(executable_file_name),
            str(start_year),
            str(year),
            str(float(lat_domain[i])),
            str(lon_domain_tmp),
            "../domain/{}".format(domain_file_name),
            "../domain/{}".format(nutrient_file_name),
            unique_job_id,
            met_data_temporary_location,
            "map",
            str(i + 1),
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
            "../domain/{}".format(domain_file_name),
            "../domain/{}".format(nutrient_file_name),
            unique_job_id,
            met_data_temporary_location,
            "map",
            str(i + 1),
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
            "EOF",
        ]
    )
    # print(run_command)
    proc = subprocess.Popen(
        [run_command], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    out, err = proc.communicate()
    return out, err