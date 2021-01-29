#!/bin/bash

# Load in config options
source ./config.sh

# Load the R module so we can use Rscript
echo "Loading R module"
module load R

# Create necessary directories and clear old data if present
mkdir -p ${SR_TIMING_DIR}/data

# Generate the list of raw timing data files
echo "Generating list of raw data files"
find ${SR_TIMING_OUTPUT_DIR} -name *_multicell.dat | sort > ${SR_TIMING_DIR}/raw_timing_data_files.txt

# Call RScript to scrape the data
echo "Running R script to scrape data"
Rscript ${SR_ROOT_DIR}/experiments/scripts/timing_scrape.R ${SR_TIMING_DIR}/raw_timing_data_files.txt ${SR_TIMING_DIR}/scraped_timing_data.csv ${SR_TIMING_SAMPLES}

echo 
