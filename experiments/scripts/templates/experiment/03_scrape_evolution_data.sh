#!/bin/bash

# Load in config options
source ./config.sh

# Load the R module so we can use Rscript
echo "Loading R module"
module load R

# Create necessary directories
SR_EVO_DATA_DIR=${SR_EVO_DIR}/data
mkdir -p ${SR_EVO_DATA_DIR}
SR_EVO_PLOT_DIR=${SR_EVO_DIR}/plots
mkdir -p ${SR_EVO_PLOT_DIR}

# Generate the list of raw evolution data files
echo "Generating list of raw data files"
find ${SR_EVO_OUTPUT_DIR} -name *__slurm.out | sort > ${SR_EVO_DATA_DIR}/raw_evolution_data_files.txt

# Call RScript to scrape the data
echo " "
echo "################"
echo "Running R script to scrape data"
Rscript ${SR_ROOT_DIR}/experiments/scripts/evolution_scrape.R ${SR_EVO_DATA_DIR}/raw_evolution_data_files.txt ${SR_EVO_DATA_DIR}/scraped_evolution_data.csv ${SR_EVO_SCRAPE_GEN_STEP} ${SR_EVO_SCRAPE_GEN_MIN} ${SR_EVO_SCRAPE_GEN_MAX} ${SR_EVO_REPS}

# Convert scraped data into .dat files that SpatialRestraint can read in
echo " "
echo "################"
echo "Running R script to convert scraped data to a more usable format"
Rscript ${SR_ROOT_DIR}/experiments/scripts/spatial_restraint_data_prep.R ${SR_EVO_DATA_DIR}/scraped_evolution_data.csv ${SR_EVO_DATA_DIR}

# Plot summary data
echo " "
echo "################"
echo "Plotting summary data"
Rscript ${SR_ROOT_DIR}/experiments/scripts/spatial_restraint_summary_plots.R ${SR_EVO_DATA_DIR}/intermediate_data_summary.csv ${SR_EVO_PLOT_DIR}

# Plot full data
echo " "
echo "################"
echo "Plotting full data"
Rscript ${SR_ROOT_DIR}/experiments/scripts/spatial_restraint_full_plots.R ${SR_EVO_DATA_DIR}/intermediate_data.csv ${SR_EVO_PLOT_DIR} ${SR_EVO_SCRAPE_GEN_MAX}
