# Scrapes *evolution* data into a usable .csv file

rm(list = ls())
library(hash)

#### BEGIN CONFIGURATION ####

args = commandArgs(trailingOnly=T)
if(length(args) != 6){
    print('Error! Must pass exactly 6 command line arguments:')
    print('1. .txt file with list of files to scrape')
    print('2. Directory to output the resulting .csv file')
    print('3. When scraping, how many generations to skip at a time')
    print('4. When scraping, what generation to start at')
    print('5. When scraping, what generation to stop at')
    print('6. Number of replicates per treatment')
    q()
}

# A .txt file, each line is a slurm output file we will scrape
file_list_filename = args[1]#'files_to_scrape.txt'
# Where to save the output .csv?
output_filename = args[2]#'data_spatial_restraint_start_75.csv'
# This one's tricky
# If we have a filename eg. /mnt/gs18/scratch/users/bob/MCSIZE_512__ONES_70/101/slurm.out
# We split it on the character /
# We need the index into that split that has the infromation (here MCSIZE_512...)
# Remember, R starts at index 1, and there is technically an empty string before /
# Therefore, for the abox example we'd need a value of 7
#filename_focal_idx = 10
# How many generations to skip between records? (1 = all records & a huge file)
gen_step = as.numeric(args[3])#100
# What are the maximum and minimum generations?
gen_min = as.numeric(args[4])#0
gen_max = as.numeric(args[5])#5000
gen_count = (gen_max - gen_min) + 1
gen_count_actual = (gen_max - gen_min) / gen_step + 1
# How many replicates exist *in each file*?
rep_count = as.numeric(args[6])#100 
cat('rep_count: ', rep_count, '\n')
cat('Gen count actual: ', gen_count_actual, '\n')

#### END CONFIGURATION ####

# Sources: 
# https://stackoverflow.com/questions/4106764/what-is-a-good-way-to-read-line-by-line-in-r

filename_vec = as.character(read.csv(file_list_filename, header = F)[,1])
# Don't initialize data until we read the first file
data = NA
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

options(warn = 1)
for(filename_idx in 1:length(filename_vec)){
    #Clean up the file via grep

    # Get the next filename, and extract all information embedded in it
    filename = filename_vec[filename_idx]
    filename_parts = strsplit(strsplit(filename, '/')[[1]][filename_focal_idx], '__')
    filename_var_hash = hash()
    filename_vars = c()
    for(filename_part in filename_parts[[1]]){
        filename_bits = strsplit(filename_part, '_')
        #if(filename_bits[[1]][1] != 'spatial'){
        if(filename_bits[[1]][1] %in% 
                c('MCSIZE', 'COST', 'GENS', 'MUT', 'POP', 'SAMPLES', 'REPS', 'ONES', 'CELLMUT', 'THRESH', 'LENGTH')){
            filename_vars = c(filename_vars, filename_bits[[1]][1])
            filename_var_hash[[filename_bits[[1]][1]]] = filename_bits[[1]][2]
        }
        if(filename_bits[[1]][1] == 'CELL'){
            filename_var_hash['CELL_MUT'] = filename_bits[[1]][3]
        }
    }
    print(filename_vars)
    print(filename_var_hash)
    cat(filename, '\n')
    if(!file.exists(filename)){
        cat('File not found! Skipping!\n')
        next    
    }
    tmp_id =  round(runif(1, 1, 120000))
    tmp_filename = paste0('tmp_file_', tmp_id,  '.out')
    tmp_filename_2 = paste0('tmp_file_', tmp_id,  '_2.out')
    system(paste0('cat ', filename, '| grep -P "^(START Treatment|#gen|\\d+,)" > ', tmp_filename))
    file_info = file.info(tmp_filename)
    if(file_info$size == 0){
        cat('File empty after grep. Skipping.\n')
        next
    }
    # Catch nasty bug where the START line is the only one that passes the previous grep
    system(paste0('cat ', filename, '| grep -P "^(#gen|\\d+,)" > ', tmp_filename_2))
    file_info_2 = file.info(tmp_filename_2)
    if(file_info_2$size == 0){
        cat('File empty after grep. Skipping.\n')
        next
    }
    # Read the file!
    fp = file(tmp_filename, open = 'r')
    rep_id = 1
    line_num = 1
    # Keep reading until we hit an empty line
    while(length(line <- readLines(fp, n = 1, warn = F)) > 0){
        # Check to see if this line is starting a new replicate
        parts_list = strsplit(line, ' ')
        if(!is.na(parts_list[[1]][1]) && parts_list[[1]][1] == 'START'){
            # Start prepping data for this replicate
            rep_id = as.numeric(parts_list[[1]][length(parts_list[[1]])])
            cat(rep_id, ' ')
            # Load *only* the data for this replicate, but do it all at once
            data_rep = read.csv(tmp_filename, skip = line_num, nrow = gen_count, header = T)
            colnames(data_rep) = c('generation', 'ave_ones', 'ave_repro_time', 'min_ones', 'max_ones', 'var_ones')
            data_rep$generation = as.numeric(data_rep$generation)
            data_rep = data_rep[data_rep$generation %% gen_step == 0,]
            data_rep$rep_id = rep_id
            for(var in filename_vars){
                data_rep[,var] = filename_var_hash[[var]]
            }
            # If data isn't initialized, initialize it here (with the sizes we now know)
            if(!data_initialized){
                data = data.frame(data = matrix(
                    nrow = rep_count * gen_count_actual * length(filename_vec), 
                    ncol = ncol(data_rep)))
                colnames(data) = colnames(data_rep)
                data_initialized = T
            }
            # Figure out where this replicate will go in the overarching data frame
            # filename_idx is 1-indexed, rep_id is 0-indexed
            start_idx = (filename_idx - 1) * rep_count * gen_count_actual + rep_id * gen_count_actual + 1
            stop_idx =  (filename_idx - 1) * rep_count * gen_count_actual + (rep_id + 1) * gen_count_actual 
            #cat(start_idx, ':', stop_idx, '\n')
            #print(data_rep$generation)
            # Insert the data!
            if(nrow(data_rep) == (stop_idx - start_idx + 1)){
                data[start_idx:stop_idx,] = data_rep
            } else {
                cat(filename, ' corrupted. Skipping the rest of it\n')
                break
            } 
            line <- readLines(fp, n = gen_count + 1, warn = F)
            line_num = line_num + 2 + gen_count # This line + header + data
        }
        # Else this is not starting a replicate, move onto the next line
        else{
            line_num = line_num + 1
        }
    }
    cat('\n')
    system(paste0('rm ', tmp_filename))
    system(paste0('rm ', tmp_filename_2))
}
# Output the data! :^)
write.csv(data, output_filename)
cat(paste0('Done! File saved to: ', output_filename, '!', '\n'))
