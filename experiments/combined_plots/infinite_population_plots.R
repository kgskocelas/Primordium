rm(list = ls())

# Load libraries
library(ggplot2)
library(cowplot)
library(khroma)

####################################
############ FINITE ################
####################################

# Load data
df_finite = read.csv('../2021_03_07__inf_population/simulation/data/size_16.csv')
df_finite = rbind(df_finite, read.csv('../2021_03_07__inf_population/simulation/data/size_32.csv'))
df_finite = rbind(df_finite, read.csv('../2021_03_07__inf_population/simulation/data/size_64.csv'))
df_finite = rbind(df_finite, read.csv('../2021_03_07__inf_population/simulation/data/size_128.csv'))
df_finite = rbind(df_finite, read.csv('../2021_03_07__inf_population/simulation/data/size_256.csv'))
df_finite = rbind(df_finite, read.csv('../2021_03_07__inf_population/simulation/data/size_512.csv'))

# Summarize data, grabbing the mean from the simulation
size_vec = c(16, 32, 64, 128, 256, 512)
df_finite_summary = data.frame(data = matrix(nrow = 0, ncol = 3))
colnames(df_finite_summary) = c('mc_size', 'generation', 'mean_ones')
for(size in size_vec){
  mask = df_finite$mc_size == size
  cat('Size: ', size, '\n\t')
  for(gen in unique(df_finite$generation)){
    df_finite_summary[nrow(df_finite_summary) + 1,] = c(size, gen, weighted.mean(0:100, df_finite[mask & df_finite$generation == gen,]$frac_of_pop))
    cat(gen, ' ')
  }
  cat('\n')
}

# Create variables to make plotting easier
df_finite_summary$restraint_value = df_finite_summary$mean_ones - 60
df_finite_summary$size_str = paste0(df_finite_summary$mc_size, 'x', df_finite_summary$mc_size)
df_finite_summary$size_factor = factor(df_finite_summary$size_str, levels = paste0(size_vec, 'x', size_vec))
text_major_size = 18
text_minor_size = 16
color_vec = as.character(khroma::color('bright')(7))
color_map = c(
  '8x8' =       '#333333',
  '16x16' =     color_vec[1],
  '32x32' =     color_vec[2],
  '64x64' =     color_vec[3],
  '128x128' =   color_vec[4],
  '256x256' =   color_vec[5],
  '512x512' =   color_vec[6],
  '1024x1024' = color_vec[7]
)

# Create directory to save plots
if(!dir.exists('./plots')){
  dir.create('./plots')
}

ggp_finite = ggplot(df_finite_summary[df_finite_summary$generation == 10000,], aes(x = size_factor, y = restraint_value ))  +
  geom_hline(aes(yintercept = 0), alpha = 0.5, linetype = 'dashed') + 
  geom_col(aes(fill = size_factor)) + 
  geom_text(aes(x = size_factor, y = restraint_value + 2, label = round(restraint_value, 1))) + 
  xlab('Organism size') + 
  ylab('Average restraint buffer') + 
  labs(color = 'Organism size') + 
  scale_fill_manual(values = color_map) + 
  scale_y_continuous(limits = c(0, 70)) +
  theme_light() + 
  theme(axis.title = element_text(size = text_major_size)) + 
  theme(axis.text = element_text(size = text_minor_size)) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 0.95, hjust = 0.9)) +
  theme(legend.title = element_text(size = text_major_size)) + 
  theme(legend.text = element_text(size = text_minor_size)) + 
  theme(panel.grid.major.x = element_blank()) + 
  theme(panel.grid.minor.x = element_blank()) + 
  theme(strip.text = element_text(size = text_minor_size, color = '#000000')) + 
  theme(strip.background = element_rect(fill = '#dddddd')) +
  theme(legend.position = 'none') 


####################################
########### INFINITE ###############
####################################

# Load data
df_infinite = read.csv('../2021_03_08__inf_genome_inf_pop/simulation/data/size_16.csv')
df_infinite = rbind(df_infinite, read.csv('../2021_03_08__inf_genome_inf_pop/simulation/data/size_32.csv'))
df_infinite = rbind(df_infinite, read.csv('../2021_03_08__inf_genome_inf_pop/simulation/data/size_64.csv'))
df_infinite = rbind(df_infinite, read.csv('../2021_03_08__inf_genome_inf_pop/simulation/data/size_128.csv'))
df_infinite = rbind(df_infinite, read.csv('../2021_03_08__inf_genome_inf_pop/simulation/data/size_256.csv'))
df_infinite = rbind(df_infinite, read.csv('../2021_03_08__inf_genome_inf_pop/simulation/data/size_512.csv'))

# Summarize data, grabbing the mean from the simulation
size_vec = c(16, 32, 64, 128, 256, 512)
df_infinite_summary = data.frame(data = matrix(nrow = 0, ncol = 3))
colnames(df_infinite_summary) = c('mc_size', 'generation', 'mean_ones')
for(size in size_vec){
  mask = df_infinite$mc_size == size
  cat('Size: ', size, '\n\t')
  for(gen in unique(df_infinite$generation)){
    df_infinite_summary[nrow(df_infinite_summary) + 1,] = c(size, gen, weighted.mean(-100:550, df_infinite[mask & df_infinite$generation == gen,]$frac_of_pop))
    cat(gen, ' ')
  }
  cat('\n')
}

# Create variables to make plotting easier
df_infinite_summary$restraint_value = df_infinite_summary$mean_ones 
df_infinite_summary$size_str = paste0(df_infinite_summary$mc_size, 'x', df_infinite_summary$mc_size)
df_infinite_summary$size_factor = factor(df_infinite_summary$size_str, levels = paste0(size_vec, 'x', size_vec))

ggp_infinite = ggplot(df_infinite_summary[df_infinite_summary$generation == 10000,], aes(x = size_factor, y = restraint_value ))  +
  geom_hline(aes(yintercept = 0), alpha = 0.5, linetype = 'dashed') + 
  geom_col(aes(fill = size_factor)) + 
  geom_text(aes(x = size_factor, y = restraint_value + 2, label = round(restraint_value, 1))) + 
  xlab('Organism size') + 
  ylab('Average restraint buffer') + 
  labs(color = 'Organism size') + 
  scale_fill_manual(values = color_map) + 
  scale_y_continuous(limits = c(0, 70)) +
  theme_light() + 
  theme(axis.title = element_text(size = text_major_size)) + 
  theme(axis.text = element_text(size = text_minor_size)) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 0.95, hjust = 0.9)) +
  theme(legend.title = element_text(size = text_major_size)) + 
  theme(legend.text = element_text(size = text_minor_size)) + 
  theme(panel.grid.major.x = element_blank()) + 
  theme(panel.grid.minor.x = element_blank()) + 
  theme(strip.text = element_text(size = text_minor_size, color = '#000000')) + 
  theme(strip.background = element_rect(fill = '#dddddd')) +
  theme(legend.position = 'none')  

####################################
########## COMBINED ################
####################################

combined_plot = plot_grid(
  ggp_finite + theme(axis.title.x = element_blank()), 
  ggp_infinite + theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank()), 
  labels = c('A', 'B'), label_size = text_major_size)
combined_plot = ggdraw(add_sub(combined_plot, "Organism Size", size = text_major_size, vpadding=grid::unit(0,"lines"),y=5, x=0.5, vjust=4.5))
plot(combined_plot)

save_plot('./plots/combined_inf_pop_plots.png',
  combined_plot,
  units = 'in', base_width = 3, base_height = 6, ncol = 2)
save_plot('./plots/combined_inf_pop_plots.pdf',
  combined_plot,
  units = 'in', base_width = 3, base_height = 6, ncol = 2)
