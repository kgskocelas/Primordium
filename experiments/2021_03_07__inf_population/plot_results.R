rm(list = ls())

# Load libraries
library(ggplot2)
library(cowplot)
library(khroma)

# Load data
df = read.csv('./simulation/data/size_16.csv')
df = rbind(df, read.csv('./simulation/data/size_32.csv'))
df = rbind(df, read.csv('./simulation/data/size_64.csv'))
df = rbind(df, read.csv('./simulation/data/size_128.csv'))
df = rbind(df, read.csv('./simulation/data/size_256.csv'))
df = rbind(df, read.csv('./simulation/data/size_512.csv'))

# Summarize data, grabbing the mean from the simulation
size_vec = c(16, 32, 64, 128, 256, 512)
df_summary = data.frame(data = matrix(nrow = 0, ncol = 3))
colnames(df_summary) = c('mc_size', 'generation', 'mean_ones')
for(size in size_vec){
  mask = df$mc_size == size
  cat('Size: ', size, '\n\t')
  for(gen in unique(df$generation)){
    df_summary[nrow(df_summary) + 1,] = c(size, gen, weighted.mean(0:100, df[mask & df$generation == gen,]$frac_of_pop))
    cat(gen, ' ')
  }
  cat('\n')
}

# Create variables to make plotting easier
df_summary$restraint_value = df_summary$mean_ones - 60
df_summary$size_str = paste0(df_summary$mc_size, 'x', df_summary$mc_size)
df_summary$size_factor = factor(df_summary$size_str, levels = paste0(size_vec, 'x', size_vec))
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
if(!dir.exists('./simulation')){
  dir.create('./simulation')
}
if(!dir.exists('./simulation/plots')){
  dir.create('./simulation/plots')
}
if(!dir.exists('./simulation/plots/final')){
  dir.create('./simulation/plots/final')
}

ggplot(df_summary[df_summary$generation <= 10000,], aes(x = generation, y = restraint_value ))  +
  geom_line(aes(color = size_factor), size = 1.0) +
  xlab('Generation') + 
  ylab('Average restraint buffer') + 
  labs(color = 'Organism size') +
  scale_color_manual(values = color_map) +
  scale_x_continuous(expand = c(0,0)) +
  theme_light() +
  theme(axis.title = element_text(size = text_major_size)) +
  theme(axis.text = element_text(size = text_minor_size)) +
  theme(legend.title = element_text(size = text_major_size)) +
  theme(legend.text = element_text(size = text_minor_size)) +
  theme(strip.text = element_text(size = text_minor_size, color = '#000000')) +
  theme(strip.background = element_rect(fill = '#dddddd')) +
  theme(legend.position = 'bottom') 
ggsave('./simulation/plots/final/simulated_density_over_time_trimmed.png', units = 'in', width = 6, height = 6) 
ggsave('./simulation/plots/final/simulated_density_over_time_trimmed.pdf', units = 'in', width = 6, height = 6)

ggplot(df_summary[df_summary$generation == 10000,], aes(x = size_factor, y = restraint_value ))  +
  geom_hline(aes(yintercept = 0), alpha = 0.5, linetype = 'dashed') + 
  geom_col(aes(fill = size_factor)) + 
  geom_text(aes(x = size_factor, y = restraint_value + 2, label = round(restraint_value, 1))) + 
  xlab('Organism size') + 
  ylab('Average restraint buffer') + 
  labs(color = 'Organism size') + 
  scale_fill_manual(values = color_map) + 
  scale_y_continuous(limits = c(0, 60)) +
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
ggsave('./simulation/plots/final/evolved_bars.png', units = 'in', width = 6, height = 6) 
ggsave('./simulation/plots/final/evolved_bars.pdf', units = 'in', width = 6, height = 6)


#df_inf_summary = read.csv('../infinite/60_40/simulation/density/combined_simulation_summary_data.csv')
#ggp_finite = ggplot(df_summary[df_summary$generation == 10000,], aes(x = size_factor, y = restraint_value ))  +
#  geom_hline(aes(yintercept = 0), alpha = 0.5, linetype = 'dashed') + 
#  geom_col(aes(fill = size_factor)) + 
#  geom_text(aes(x = size_factor, y = restraint_value + 2, label = round(restraint_value, 1))) + 
#  xlab('Organism size') + 
#  ylab('Average restraint buffer') + 
#  labs(color = 'Organism size') + 
#  #scale_color_manual(values = color_map) + 
#  scale_y_continuous(limits = c(0, 60)) +
#  theme_light() + 
#  theme(axis.title = element_text(size = text_major_size)) + 
#  theme(axis.text = element_text(size = text_minor_size)) + 
#  theme(axis.text.x = element_text(angle = 60, vjust = 0.95, hjust = 0.9)) +
#  theme(legend.title = element_text(size = text_major_size)) + 
#  theme(legend.text = element_text(size = text_minor_size)) + 
#  theme(panel.grid.major.x = element_blank()) + 
#  theme(panel.grid.minor.x = element_blank()) + 
#  theme(strip.text = element_text(size = text_minor_size, color = '#000000')) + 
#  theme(strip.background = element_rect(fill = '#dddddd')) +
#  theme(legend.position = 'none') #
#  #ggsave('./density_3/test.png', units = 'in', width = 3, height = 6) +
#  #ggsave('./density_3/test.pdf', units = 'in', width = 3, height = 6)
#
#df_inf_summary$restraint_value = df_inf_summary$mean_ones 
#df_inf_summary$size_str = paste0(df_inf_summary$mc_size, 'x', df_inf_summary$mc_size)
#size_vec = c(16, 32, 64, 128, 256, 512, 1024)
#df_inf_summary$size_factor = factor(df_inf_summary$size_str, levels = paste0(size_vec, 'x', size_vec))
#ggp_inf = ggplot(df_inf_summary[df_inf_summary$generation == 10000 & df_inf_summary$mc_size != 1024,], aes(x = size_factor, y = restraint_value ))  +
#  geom_hline(aes(yintercept = 0), alpha = 0.5, linetype = 'dashed') + 
#  geom_col(aes(fill = size_factor)) + 
#  geom_text(aes(x = size_factor, y = restraint_value + 2, label = round(restraint_value, 1))) + 
#  #geom_text(aes(x = size_factor, y = restraint_value + 1, label = round(restraint_value, 2))) + 
#  xlab('Organism size') + 
#  ylab('Average restraint buffer') + 
#  labs(color = 'Organism size') + 
#  #scale_color_manual(values = color_map) + 
#  scale_y_continuous(limits = c(0, 60)) +
#  theme_light() + 
#  theme(axis.title = element_text(size = text_major_size)) + 
#  theme(axis.text = element_text(size = text_minor_size)) + 
#  theme(axis.text.x = element_text(angle = 60, vjust = 0.95, hjust = 0.9)) + 
#  theme(legend.title = element_text(size = text_major_size)) + 
#  theme(legend.text = element_text(size = text_minor_size)) + 
#  theme(panel.grid.major.x = element_blank()) + 
#  theme(panel.grid.minor.x = element_blank()) + 
#  theme(strip.text = element_text(size = text_minor_size, color = '#000000')) + 
#  theme(strip.background = element_rect(fill = '#dddddd')) +
#  theme(legend.position = 'none')
#  #ggsave('./density_3/test.png', units = 'in', width = 3, height = 6) +
#  #ggsave('./density_3/test.pdf', units = 'in', width = 3, height = 6)
#
#combined_plot = plot_grid(ggp_finite, ggp_inf, labels = c('A', 'B'), label_size = 18)
#plot(combined_plot)
#
#save_plot('./density_3/simulated_combined_test.png',
#  combined_plot,
#  units = 'in', base_width = 3, base_height = 6, ncol = 2)
#save_plot('./density_3/simulated_combined_test.pdf',
#  combined_plot,
#  units = 'in', base_width = 3, base_height = 6, ncol = 2)

