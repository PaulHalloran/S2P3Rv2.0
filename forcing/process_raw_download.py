import subprocess
import os
import glob
import datetime
import pandas as pd
from os.path import exists
import logging

base_dir = '/massive/ph290/cmip6/'
model_directory_name = 'UKESM1'
model_name='UKESM1-0-LL'
variables=['tas','uas','vas','psl','hurs','rsds','rlds']

dict={}
dict[model_name]={}
dict[model_name]['experiments'] = ['ssp585']
dict[model_name]['grid_labels']=['gn']
dict[model_name]['variant_labels']=['r1i1p1f2']

Log_Format = "%(levelname)s %(asctime)s - %(message)s"

logging.basicConfig(filename = "process_raw_download.log",
                    filemode = "w",
                    format = Log_Format,
                    level = logging.ERROR)

logger = logging.getLogger()
logger.error('##################################################################')
logger.error('########  New run                                          #######')
logger.error('##################################################################')

for model in dict.keys():
  for experiment in dict[model]['experiments']:
    for grid_label in dict[model]['grid_labels']:
      for variant_label in dict[model]['variant_labels']:
        for var in variables:
            print('Attempting to process: '+var+'_day_'+model+'_'+experiment+'_'+variant_label+'_'+grid_label+'.nc')
            logger.error('Attempting to process: '+var+'_day_'+model+'_'+experiment+'_'+variant_label+'_'+grid_label+'.nc')
            if not exists(base_dir+'/'+experiment+'/'+model_directory_name+'/'+var+'_day_'+model+'_'+experiment+'_'+variant_label+'_'+grid_label+'.nc'):
                ## Check that all of the required files are present
                years_list = []
                files = glob.glob(base_dir+'/'+experiment+'/'+model_directory_name+'/'+var+'_day_'+model+'_'+experiment+'_'+variant_label+'_'+grid_label+'_????????-????????.nc')
                for file in files:
                    years = file.split('/')[-1].split('_')[-1].split('.')[0].split('-')
                    years_list.append([datetime.datetime.strptime(year,'%Y%m%d') for year in years])

                ### Sort the list (based on 1st year in each subarray)
                years_list = sorted(years_list, key=lambda x: x[0])

                ### Identify and display the min and max year to check that teh required range have been downloaded
                years_flat =  [item for sublist in years_list for item in sublist]
                first_year = min(years_flat)
                last_year = max(years_flat)
                print('Year range in data is: ',first_year,last_year)
                logger.error('Year range in data is: ',first_year,last_year)

                ### Check that all files are present, i.e. that there are no year ranges that have not downloaded
                for i,dummy in enumerate(years_list[1::]):
                    date_diff_between_sequential_files = years_list[i+1][0] - years_list[i][-1]
                    if date_diff_between_sequential_files > datetime.timedelta(days = 20):
                        print('There appear to be missing files. Check dates in filelist: ')
                        print(files)
                        logger.error('There appear to be missing files. Check dates in filelist: ')
                        logger.error(files)
                        exit()


                print('necessary files are present')
                print('merging files')

                ## Files are pressent, so merge into single file per variable
                subprocess.call(['cdo mergetime '+base_dir+'/'+experiment+'/'+model_directory_name+'/'+var+'_day_'+model+'_'+experiment+'_'+variant_label+'_'+grid_label+'_????????-????????.nc  '+base_dir+'/'+experiment+'/'+model_directory_name+'/'+var+'_day_'+model+'_'+experiment+'_'+variant_label+'_'+grid_label+'.nc'], shell=True)
                print('removing unmerged files')
                files_to_remove = glob.glob(base_dir+'/'+experiment+'/'+model_directory_name+'/'+var+'_day_'+model+'_'+experiment+'_'+variant_label+'_'+grid_label+'_????????-????????.nc')
                # for file_to_remove in files_to_remove:
                #     os.remove(file_to_remove)
            else:
                print('Output files already exists: '+var+'_day_'+model+'_'+experiment+'_'+variant_label+'_'+grid_label+'.nc')
