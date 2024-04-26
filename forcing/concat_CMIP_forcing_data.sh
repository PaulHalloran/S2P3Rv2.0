#!/bin/bash -l

###############################################################################
# concat_CMIP_forcing_data.sh
# by Joshua Wiggs <joshua.wiggs@metoffice.gov.uk>
###############################################################################

usage(){
  cat << EOF
  Usage: $0 [-e] [-m] [-s] [-t]

  Application for creating the CMIP6 forcing data required to produce 
  meteorological data to drive the S2P3 V2.0 Shelf Seas model.

  Options:
  -e| Select the CMIP6 experiment 
  -m| Select the CMIP6 model
  -s| Select the CMIP6 ensemble member
  -t| Set the top level directory for forcing data

EOF
}

#Add flag functionality -e -m -s -t
while getopts "e:m:s:t:" o; do
  case ${o} in
    e)
      EXPERIMENT=${OPTARG}
      ;;
    m)
      MODEL=${OPTARG}
      ;;
    s)
      ENSEMBLE_MEMBER=${OPTARG}
      ;;
    t)
      TARGET_DIRECTORY=${OPTARG}
      ;;
    ?)
      usage
      exit 0
      ;;
  esac
done
shift $((OPTIND-1))

# Set which CMIP6 variables to produce forcing data from
VARIABLES=("hurs" "psl" "rlds" "rsds" "tas" "uas" "vas" "sfcWind")

# Activate conda environment with CDO installed
conda activate /home/h06/jwiggs/.conda/envs/cdo_env

# Ensure correct folders are present in target directory to store data
OUTPUT_DIRECTORY="${TARGET_DIRECTORY}/${EXPERIMENT}/${MODEL}"
mkdir -p $OUTPUT_DIRECTORY

# Loop through each variable to produce forcing data
for ITER in ${!VARIABLES[@]}; do
  VARIABLE=${VARIABLES[$ITER]}
  echo "Starting ${VARIABLE}..."

  # Get list of CMIP6 files using managecmip
  CMIP_FILES=($(managecmip list-local -M CMIP6 -t day -m $MODEL \
    -e $EXPERIMENT -s $ENSEMBLE_MEMBER -v $VARIABLE))
  HISTORIC_CMIP_FILES=($(managecmip list-local -M CMIP6 -t day -m $MODEL \
    -e historical -s $ENSEMBLE_MEMBER -v $VARIABLE))

  if [[ ${#CMIP_FILES[@]} == 0 ]]; then
    echo "No CMIP files found for experiment ${EXPERIMENT}"
    continue
  fi

  if [[ ${#HISTORIC_CMIP_FILES[@]} == 0 ]]; then
    echo "No historical CMIP files found for ${MODEL}"
    continue
  fi
  
  # Pharse historical filenames to select one starting in 1950
  for FILE_ITER in ${!HISTORIC_CMIP_FILES[@]}; do
    FILE=${HISTORIC_CMIP_FILES[FILE_ITER]}
    IFS="/" read -ra SPLIT_FILE <<< "$FILE"
    FILENAME=${SPLIT_FILE[-1]}
    IFS="_" read -ra SPLIT_FILENAME <<< "$FILENAME"
    FILE_YEAR=${SPLIT_FILENAME[-1]}
    if [[ "${FILE_YEAR:0:17}" == "19500101-20141230" ]]; then
        HISTORICAL_CMIP_FILE=$FILE
    fi
  done

  # Construct filename for merged CMIP6 forcing data
  IFS="/" read -ra FNAME <<< "${CMIP_FILES[0]}"
  FILEID=$(echo "${FNAME[-1]}" | awk -F '_' '{print $2"_"$3"_"$4"_"$5}')
  NEW_FNAME="${OUTPUT_DIRECTORY}/${VARIABLE}_${FILEID}_all.nc"

  # Take CMIP6 files and merge into single forcing dataset
  cdo mergetime "$HISTORICAL_CMIP_FILE" "${CMIP_FILES[@]}" "$NEW_FNAME"

done
