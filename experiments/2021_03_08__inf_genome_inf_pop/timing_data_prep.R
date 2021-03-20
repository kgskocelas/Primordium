rm(list = ls())

# Load libraries
library(ggplot2)
library(ggridges)
library(dplyr)

# Load the data and trim
df = read.csv('./timing/data/scraped_timing_data.csv')
df = df[df$cells_side != 8,]

# Group and summarize data
df_grouped = dplyr::group_by(df, cells_side, ancestor_1s)
df_summary = dplyr::summarize(df_grouped, mean_time = mean(rep_time), n = n(), sd_time = sd(rep_time))
# Make sure no data is missing
df_summary[df_summary$n != 100,]

# Normalize data by lowest restraint level
df_summary$mean_time_norm = 0
df_summary$sd_time_norm = 0
for(cells_side in unique(df_summary$cells_side)){
  mask = df_summary$cells_side == cells_side
  norm_factor = df_summary[mask & df_summary$ancestor_1s == -100,]$mean_time
  df_summary[mask,]$mean_time_norm = df_summary[mask,]$mean_time / norm_factor
  df_summary[mask,]$sd_time_norm = df_summary[mask,]$sd_time / norm_factor
}

# Calculate the reproduction rate, simply 1 / normalized reproduction time
df_summary$mean_rep_rate = 1 / df_summary$mean_time_norm

# Calculate the fitness advantage of this restraint buffer value compared to one less bit of restraint
df_summary$rep_rate_diff = 0
for(cells_side in unique(df_summary$cells_side)){
  mask = df_summary$cells_side == cells_side
  rep_rates = df_summary[mask,]$mean_rep_rate
  diffs = rep_rates[2:length(rep_rates)] - rep_rates[1:(length(rep_rates)-1)]
  df_summary[mask,]$rep_rate_diff = c(0, diffs)
}

# Save data
write.csv(df_summary, 'timing/data/timing_data_summary.csv')
