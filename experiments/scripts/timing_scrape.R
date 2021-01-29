# Scrapes multicell reproduction data into a usable .csv file
# NOTE: This is *NOT* for scraping it to be loaded into the Spatial Restraint app
    # If that's what you need, run this, and then run convert_timing_to_dat.R

rm(list = ls())
library(hash)

#### BEGIN CONFIGURATION ####

args = commandArgs(trailingOnly=T)
if(length(args) != 3){
    print('Error! Must pass exactly 3 command line arguments:')
    print('File list, output path, and the number of replicates per file')
    q()
}

# A .txt file, each line is a multicell.dat file we will scrape
file_list_filename = args[1]#'file_lists/512_timing_lower.txt'
# Where to save the output .csv?
output_filename = args[2]#'data_timing_512_lower.csv'
# This one's tricky
# If we have a filename eg. /mnt/gs18/scratch/users/bob/MCSIZE_512__ONES_70/101/multicell.dat
# We split it on the character /
# We need the index into that split that has the infromation (here MCSIZE_512...)
# Remember, R starts at index 1, and there is technically an empty string before /
# Therefore, for the abox example we'd need a value of 7
# filename_focal_idx = as.numeric(args[3])
# The number of multicell timing runs (replicates?) in each file
runs_per_file = as.numeric(args[3])

#### END CONFIGURATION ####

## Load vector of filenames, prep main data frame
filename_vec = as.character(read.csv(file_list_filename, header = F)[,1])
data = data.frame(data = matrix(nrow = length(filename_vec) * runs_per_file, ncol = 10))
data_initialized = F

# Find the focal index for filenames -- the directory that containts metadata
filename_focal_idx = NA
example_filename = filename_vec[1]
filename_part_vec = strsplit(example_filename, '/')[[1]]
for(filename_part_idx in 1:length(filename_part_vec)){
    filename_part = filename_part_vec[filename_part_idx]
    if(grepl('MCSIZE_', filename_part, fixed = T)){
        filename_focal_idx = filename_part_idx
    }
}
if(is.na(filename_focal_idx)){
    print('Error! Unable to find filename focal index in filename:')
    print(example_filename)
    q()
}
print(paste0('Assuming focal index is: ', filename_focal_idx))

print('Printing filename index as it is scraped...')
print(paste0('Out of: ', length(filename_vec)))
options(warn = 1)
# Go through each file specified in the list
for(filename_idx in 1:length(filename_vec)){
    filename = filename_vec[filename_idx]
    # Extract as much data as possible about the data
    filename_parts = strsplit(strsplit(filename, '/')[[1]][filename_focal_idx], '__')
    filename_var_hash = hash()
    filename_vars = c()
    for(filename_part in filename_parts[[1]]){
        filename_bits = strsplit(filename_part, '_')
        if(filename_bits[[1]][1] != 'spatial'){
            filename_vars = c(filename_vars, filename_bits[[1]][1])
            filename_var_hash[[filename_bits[[1]][1]]] = filename_bits[[1]][2]
        }
    }
    #cat(filename, ' ', filename_idx, '\n')
    cat(filename_idx, ' ')
    # Read in the current file
    data_file = read.csv(filename)
    # If we haven't set the main data frame's columns do that now.
    if(!data_initialized){
        colnames(data) = c(colnames(data_file)[1:9], 'rep_time')
        data_initialized = T
    }
    # Different runs come as different columns, but we stick them in the data as different rows
    start_idx = (filename_idx - 1) * runs_per_file + 1
    stop_idx = (filename_idx) * runs_per_file 
    data[start_idx:stop_idx, 1:9] = data_file[1,1:9]
    vals = data_file[1,10:(10+runs_per_file-1)]
    data[start_idx:stop_idx, 10] = as.numeric(vals)
    #print(data[start_idx:stop_idx, 10])
    #for(run_idx in 1:runs_per_file){
    #    data[(filename_idx - 1) * runs_per_file + run_idx,] = 
    #        c(data_file[1,1:9], data_file[1,9 + run_idx])
    #}
}
cat('\n')
# Write the data to disk! :^)
write.csv(data, output_filename)
cat(paste0('Done! File saved to: ', output_filename, '!', '\n'))
