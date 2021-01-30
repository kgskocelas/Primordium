rm(list = ls())
library(dplyr)

#### BEGIN CONFIGURATION ####

args = commandArgs(trailingOnly=T)
if(length(args) != 2){
    print('Error! Must pass exactly 6 command line arguments:')
    print('1. Filename of input (scraped .csv data file)')
    print('2. Output directory')
    q()
}
# Path to the csv data file
filename_data = args[1]#'./raw_data/data_spatial_restraint_finished.csv'
# Where to save cleaned data
output_dir = args[2]
#### END CONFIGURATION ####

#### Load the data
print('Loading data, this may take a moment.')
data = read.csv(filename_data)
print('Done loading data!')
data$NUM = NA

# Note: Sometimes you may need to load in multiple data files and merge them. 
# Here's a quick example:
#data_2 = read.csv('data_spatial_restraint_start_50_512_attempt_2.csv')
#data_2$ONES = 50
#data_2$NUM = NA
#data = rbind(data, data_2)


#### Clean up the data just a little
# Remove lines that all NA
data = data[!is.na(data$MCSIZE),]
# Truncate data to first 5k generations
#data = data[as.numeric(data$generation) <= 5000,]
if(!('ONES' %in% colnames(data))){
  print('Error!')
  print('Data should have column named ONES that containts the number of ones in the ancestor')
  print('Fix: data$ONES = 50, where data is your data.frame and 50 is the correct value')
  quit()
}
if(!('NUM' %in% colnames(data))){
  print('Error!')
  print('Data should have column named NUM, though the value is not super important')
  print('(It\'s needed for merging data frames')
  print('Fix: data$NUM = NA, where data is your data.frame')
  quit()
}

# Save off the cleaned data
filename_clean = paste0(output_dir, '/intermediate_data.csv')
write.csv(data, filename_clean)
print(paste0('Cleaned data saved to ', filename_clean, '!'))

#### Summary data prep
# Group data my condition
data_grouped = dplyr::group_by(data, generation, MCSIZE, COST, ONES)
# Take summary stats for each condition
data_summary = dplyr::summarise(data_grouped, grand_mean_ones = mean(ave_ones), sd_ave_ones = sd(ave_ones), n = n())

# 95% confidence interval
z = qnorm(0.95)
data_summary$offset_95 = z * (data_summary$sd_ave_ones / sqrt(data_summary$n)) 

# Save off summary data
filename_summary = paste0(output_dir, '/intermediate_data_summary.csv')
write.csv(data_summary, filename_summary)
print(paste0('Cleaned data saved to ', filename_summary, '!'))


print('Done!')
