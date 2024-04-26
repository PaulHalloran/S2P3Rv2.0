# For cases where u and v are on the staggered grid
# regrid u and v and store all times as a single file
import os
import iris
from iris.util import equalise_attributes

rootpath = '/project/champ/data/CMIP6/ScenarioMIP/MOHC/UKESM1-0-LL/ssp245/r1i1p1f2/day'
vars = ['uas', 'vas']
target_path='/project/ciid/projects/WISERAP_Sunda_Shelf/ssp245/UKESM1-0-LL'

tasfile = os.path.join(target_path,'tas*.nc')
tas = iris.load_cube(tasfile)

# loop over files, regrid and save to target_dir
for var in vars:
	infiles = os.path.join(rootpath, var, 'gn/v20190715', f'{var}*.nc')
	print(infiles)
	wind = iris.load(infiles)
	# remove attributes which don't match
	equalise_attributes(wind)
	wind=wind.concatenate_cube()
	wind= wind.regrid(tas,iris.analysis.Linear())
	outfile = os.path.join(target_path, f'{var}_day_UKESM1-0-LL_ssp245_r1i1p1f2_all.nc')
	iris.save(wind,outfile)

