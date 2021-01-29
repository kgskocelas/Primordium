# Loads in data output by timing_scrape.R
# Outputs it into .dat files readable by the SpatialRestraint app

rm(list = ls())
library(dplyr)


#### BEGIN CONFIGURATION ####

args = commandArgs(trailingOnly=T)
if(length(args) != 2){
    print('Error! Must pass exactly 2 command line arguments:')
    print('Input file of scraped data and the root directory for output')
    q()
}

# Where is the .csv output by timing_scrape.R? 
input_filename = args[1]#'~/research/rogue_cell/SpatialRestraint/scraped_data/512_timing_all.csv'
# Where to save output data? 
    # NOTE: this is formatted such that data will be saved to subdirectories of the path specified
        # For example, if output_dir = '/foo/bar/', the script may save it in /for/bar/512/50.dat
    # NOTE: this assumes those subdirectories ALREADY exist, else it will fail
output_dir= args[2]#'~/research/rogue_cell/SpatialRestraint/distribution_samples/'

#### END CONFIGURATION ####


# Load in data
data = read.csv(input_filename)
data = data[!is.na(data$cells_side),]

# Iterate through each size of multicell (denoted by the number of cells on one side)
for(cells_side in sort(unique(data$cells_side))){
    cat('Formatting data for' , cells_side, 'cells per side.\n')
    if(!dir.exists(file.path(output_dir, cells_side))){
        dir.create(file.path(output_dir, cells_side))
    }
    data_size = data[data$cells_side == cells_side,]
    # Within that size of multicell, iterate through all counts of ones 
    cat('Num ones: ')
    for(num_ones in unique(data_size$ancestor_1s)){
        cat(num_ones, ' ')
        # Save that data!
        data_ones = data_size[data_size$ancestor_1s == num_ones,]
        filename = paste0(output_dir, '/', cells_side, '/', num_ones, '.dat')
        write(data_ones$rep_time, filename, ncolumns = 1)
    }
    cat('\n')
}

