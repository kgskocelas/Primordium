rm(list = ls())

# Load libraries
library(ggplot2)
library(ggridges)
library(dplyr)

num_samples = 3000
num_batches = 30 
batch_size = num_samples / num_batches
reference_restraint_value = -100

# Load the data and trim
df = read.csv('./timing/data/scraped_timing_data_16_256.csv')
df = rbind(df, read.csv('./timing/data/scraped_timing_data_512.csv'))
df = df[df$cells_side != 8,]
df = df[!is.na(df$rep_time),]
df$id = 0
df$batch_id = 0

for(cells_side in unique(df$cells_side)){
  cat(paste0('MCSIZE: ', cells_side))
  cells_side_mask = df$cells_side == cells_side
  for(ancestor_1s in unique(df$ancestor_1s)){
    cat(paste0(' ', ancestor_1s))
    df[cells_side_mask & df$ancestor_1s == ancestor_1s,]$id = 1:nrow(df[cells_side_mask & df$ancestor_1s == ancestor_1s,])
    #for(batch_id in 1:num_batches){
    #  df[cells_side_mask & df$ancestor_1s == ancestor_1s,][((batch_id - 1) * batch_size + 1):((batch_id) * batch_size),]$batch_id = batch_id
    #}
  }
  cat('\n')
}

df$batch_id = floor((df$id - 1) / batch_size) + 1


# Group and summarize data
df_grouped = dplyr::group_by(df, cells_side, ancestor_1s, batch_id)
df_summary = dplyr::summarize(df_grouped, mean_time = mean(rep_time), n = n(), sd_time = sd(rep_time))
# Make sure no data is missing
df_summary[df_summary$n != batch_size,]

# Normalize data by lowest restraint level
df_summary$mean_time_norm = 0
df_summary$sd_time_norm = 0
for(cells_side in unique(df_summary$cells_side)){
  for(batch_id in 1:num_batches){
    mask = df_summary$cells_side == cells_side & df_summary$batch_id == batch_id
    norm_factor = df_summary[mask & df_summary$ancestor_1s == -reference_restraint_value,]$mean_time
    df_summary[mask,]$mean_time_norm = df_summary[mask,]$mean_time / norm_factor
    df_summary[mask,]$sd_time_norm = df_summary[mask,]$sd_time / norm_factor
  }
}

# Calculate the reproduction rate, simply 1 / normalized reproduction time
df_summary$mean_rep_rate = 1 / df_summary$mean_time_norm

# Calculate the fitness advantage of this restraint buffer value compared to one less bit of restraint
df_summary$rep_rate_diff = 0
for(cells_side in unique(df_summary$cells_side)){
  for(batch_id in 1:num_batches){
    mask = df_summary$cells_side == cells_side & df_summary$batch_id == batch_id
    rep_rates = df_summary[mask,]$mean_rep_rate
    diffs = rep_rates[2:length(rep_rates)] - rep_rates[1:(length(rep_rates)-1)]
    df_summary[mask,]$rep_rate_diff = c(0, diffs)
  }
}

# Save data
write.csv(df_summary, 'timing/data/timing_data_summary.csv')

