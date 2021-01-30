#!/bin/bash

# Load in config options
source ./config.sh

# Load the R module so we can use Rscript
echo "Loading R module"
module load R


# Generate the list of raw timing data files
echo "Generating list of raw data files"
find ${SR_EVO_OUTPUT_DIR} -name *__slurm.out | sort > ${SR_EVO_DIR}/raw_evolution_data_files.txt

# Call RScript to scrape the data
echo "Running R script to scrape data"
Rscript ${SR_ROOT_DIR}/experiments/scripts/evolution_scrape.R ${SR_EVO_DIR}/raw_evolution_data_files.txt ${SR_EVO_DIR}/scraped_evoluion_data.csv ${SR_EVO_SCRAPE_GEN_STEP} ${SR_EVO_SCRAPE_GEN_MIN} ${SR_EVO_SCRAPE_GEN_MAX} ${SR_EVO_REPS}

## Create necessary directories
#mkdir -p ${SR_EVO_DIR}/data
#
## Convert scraped data into .dat files that SpatialRestraint can read in
#echo "Running R script to convert scraped data to SR usable .dat files"
#Rscript ${SR_ROOT_DIR}/experiments/scripts/convert_timing_to_dat.R ${SR_EVO_DIR}/scraped_timing_data.csv ${SR_EVO_DIR}/data
