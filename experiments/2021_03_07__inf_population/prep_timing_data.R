rm(list = ls())

# Load libraries
library(ggplot2)
library(ggridges)
library(dplyr)

# Load in data, if recalc_summary is false we just load the summary data from file
recalc_summary = T
if(recalc_summary){
  df = read.csv('./timing/data/scraped_timing_data.csv')
  df = df[df$cells_side != 8 & df$cells_side != 1024,]
  # Group and summarize the data, we mostly just use the mean times  
  df_grouped = dplyr::group_by(df, cells_side, ancestor_1s)
  df_summary = dplyr::summarize(df_grouped, mean_time = mean(rep_time), n = n(), sd_time = sd(rep_time))
  df_summary[df_summary$n != 10000,]
}else{
  df_summary = read.csv('timing/data/timing_data_summary.csv')
}
# Add variables to make plotting easier
df_summary$size_str = paste0(df_summary$cells_side, 'x', df_summary$cells_side)
df_summary$size_factor = factor(df_summary$size_str, levels = c('8x8', '16x16', '32x32', '64x64', '128x128', '256x256', '512x512'))
df_summary$restraint_value = df_summary$ancestor_1s - 60
text_major_size = 18
text_minor_size = 16 
color_vec = as.character(khroma::color('bright')(7))
color_map = c(
  '16x16' =     color_vec[1],
  '32x32' =     color_vec[2],
  '64x64' =     color_vec[3],
  '128x128' =   color_vec[4],
  '256x256' =   color_vec[5],
  '512x512' =   color_vec[6],
  '1024x1024' = color_vec[7] 
)

# Create directory to save plots
if(!dir.exists('./timing/plots')){
  dir.create('./timing/plots')
}

# Calculate the normalized replication time
# Divide through by replication time at restraint buffer zero
df_summary$mean_time_norm = 0
df_summary$sd_time_norm = 0
for(cells_side in unique(df_summary$cells_side)){
  mask = df_summary$cells_side == cells_side
  norm_factor = df_summary[mask & df_summary$ancestor_1s == 0,]$mean_time
  df_summary[mask,]$mean_time_norm = df_summary[mask,]$mean_time / norm_factor
  df_summary[mask,]$sd_time_norm = df_summary[mask,]$sd_time / norm_factor
}

# To turn replication rate into a fitness score, we invert replication rate
df_summary$mean_rep_rate = 1 / df_summary$mean_time_norm

# Calculate the _difference_ in fitness
# At n ones, this is fit(n) - fit(n - 1)
df_summary$rep_rate_diff = 0
for(cells_side in unique(df_summary$cells_side)){
  mask = df_summary$cells_side == cells_side
  rep_rates = df_summary[mask,]$mean_rep_rate
  diffs = rep_rates[2:length(rep_rates)] - rep_rates[1:(length(rep_rates)-1)]
  df_summary[mask,]$rep_rate_diff = c(0, diffs)
}

# Save off the summary data to use in the simulation
write.csv(df_summary, 'timing/data/timing_data_summary.csv')

# Plot the replication rates of each organism size
  # x-axis = restraint buffer
  # y-axis = replication rate
  # color = organism size
ggplot(df_summary, aes(x = restraint_value, y = mean_rep_rate, color = as.factor(size_factor))) + 
  geom_vline(aes(xintercept = 0), alpha = 0.5, linetype = 'dashed') +
  geom_line(size = 1.0) + 
  xlab('Restraint buffer') + 
  ylab('Mean fitness') + 
  labs(color = 'Organism size') +
  scale_color_manual(values = color_map) +
  theme_light() +
  theme(panel.grid.major.x = element_blank()) +
  theme(axis.title = element_text(size = text_major_size)) +
  theme(axis.text = element_text(size = text_minor_size)) +
  theme(legend.title = element_text(size = text_major_size)) +
  theme(legend.text = element_text(size = text_minor_size)) +
  theme(legend.position = 'bottom') 
ggsave('timing/plots/replication_rate.png', units = 'in', width = 6, height = 6) 
ggsave('timing/plots/replication_rate.pdf', units = 'in', width = 6, height = 6)

# Plot the replication rates of each organism size
  # x-axis = restraint buffer
  # y-axis = change in fitness
  # color = organism size
ggplot(df_summary, aes(x = ancestor_1s - 60, y = rep_rate_diff, color = as.factor(size_factor))) +
  geom_vline(aes(xintercept = 0), alpha = 0.5, linetype = 'dashed') +
  geom_line(size = 1.0) + 
  xlab('Simulated restraint buffer') + 
  ylab('Change in fitness') +
  labs(color = 'Organism size') +
  scale_color_manual(values = color_map) +
  theme_light() +
  theme(panel.grid.major.x = element_blank()) +
  theme(axis.title = element_text(size = text_major_size)) +
  theme(axis.text = element_text(size = text_minor_size)) +
  theme(legend.title = element_text(size = text_major_size)) +
  theme(legend.text = element_text(size = text_minor_size)) +
  theme(legend.position = 'bottom') 
ggsave('timing/plots/replication_rate_diff.png', units = 'in', width = 6, height = 6) 
ggsave('timing/plots/replication_rate_diff.pdf', units = 'in', width = 6, height = 6)

# Plot the replication rates of each organism size
  # x-axis = restraint buffer
  # y-axis = change in fitness
  # color/subplot = organism size
ggplot(df_summary[df_summary$ancestor_1s >= 40,], aes(x = ancestor_1s - 60, y = rep_rate_diff, color = as.factor(size_factor))) +
  geom_hline(aes(yintercept = 0), alpha = 0.5, linetype = 'dashed') +
  geom_vline(aes(xintercept = 0), alpha = 0.5, linetype = 'dashed') +
  geom_line(size = 1.0) + 
  xlab('Simulated restraint bufffer') + 
  ylab('Change in fitness') +
  labs(color = 'Organism\n  size') +
  scale_color_manual(values = color_map) +
  theme_light() +
  facet_wrap(vars(size_factor), ncol = 2) +
  theme(panel.grid.major.x = element_blank()) +
  theme(axis.title = element_text(size = text_major_size)) +
  theme(axis.text = element_text(size = text_minor_size)) +
  theme(legend.title = element_text(size = text_major_size)) +
  theme(legend.text = element_text(size = text_minor_size)) +
  theme(strip.text = element_text(size = text_minor_size, color = '#000000')) +
  theme(strip.background = element_rect(fill = '#dddddd')) +
  theme(legend.position = 'none') 
ggsave('timing/plots/replication_rate_diff_trimmed_two_col.png', units = 'in', width = 6, height = 6) 
ggsave('timing/plots/replication_rate_diff_trimmed_two_col.pdf', units = 'in', width = 6, height = 6)

