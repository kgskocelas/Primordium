rm(list = ls())

# Load relevant libraries
library(dplyr)
library(ggplot2)
library(ggridges)
library(scales)
library(khroma)

# Load the data
df = read.csv('./evolution/data/scraped_evolution_data_length_50.csv')
df = rbind(df, read.csv('./evolution/data/scraped_evolution_data_length_200.csv'))
df = rbind(df, read.csv('./evolution/data/scraped_evolution_data_length_100.csv'))
df = rbind(df, read.csv('./evolution/data/scraped_evolution_data_length_400.csv'))
df = rbind(df, read.csv('./evolution/data/scraped_evolution_data_length_25.csv'))

# Create directory to save plots
if(!dir.exists('./evolution/plots')){
  dir.create('./evolution/plots')
}

# Trim off NAs (artifiacts of how we scraped the data) and trim to only have gen 10,000
df2 = df[!is.na(df$MCSIZE) & df$generation == 10000,]
# Ignore data for size 8x8 and 1024x1024
df2 = df2[df2$MCSIZE != 8 & df2$MCSIZE != 1024,]

# Group the data by size and summarize
data_grouped = dplyr::group_by(df2, MCSIZE, LENGTH)
data_summary = dplyr::summarize(data_grouped, mean_ones = mean(ave_ones), n = dplyr::n())

# Print any data that is missing!
if(sum(data_summary$n != 100) != 0){
  cat('Error! Missing data detected!\n')
  print(data_summary[data_summary$n != 100,])
}

## Set variables to make plotting easier
# Calculate restraint value (x - 60% of the genome length)
df2$restraint_value = df2$ave_ones - (df2$LENGTH * 0.6)
# Make a nice, clean factor for size
df2$size_str = paste0(df2$MCSIZE, 'x', df2$MCSIZE)
df2$size_factor = factor(df2$size_str, levels = c('16x16', '32x32', '64x64', '128x128', '256x256', '512x512', '1024x1024'))
df2$size_factor_reversed = factor(df2$size_str, levels = rev(c('16x16', '32x32', '64x64', '128x128', '256x256', '512x512', '1024x1024')))
df2$length_str = paste0(df2$LENGTH, '-bit')
df2$length_factor = factor(df2$length_str, levels = c('25-bit', '50-bit', '100-bit', '200-bit', '400-bit'))
data_summary$size_str = paste0(data_summary$MCSIZE, 'x', data_summary$MCSIZE)
data_summary$size_factor = factor(data_summary$size_str, levels = c('16x16', '32x32', '64x64', '128x128', '256x256', '512x512', '1024x1024'))
data_summary$length_str = paste0(data_summary$LENGTH, '-bit')
data_summary$length_factor = factor(data_summary$length_str, levels = c('25-bit', '50-bit', '100-bit', '200-bit', '400-bit'))
# Create a map of colors we'll use to plot the different organism sizes
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
# Set the sizes for text in plots
text_major_size = 18
text_minor_size = 16 

# Plot the number of replicates for each organism size
ggplot(data_summary, aes(x = size_factor, y = n)) +
  geom_col(aes(fill = size_factor)) +
  geom_text(aes(y = n + 5, label = n)) +
  scale_fill_manual(values = color_map) +
  theme(legend.position = 'none') +
  theme(panel.grid.major.x = element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(axis.title = element_text(size = text_major_size)) +
  theme(axis.text = element_text(size = text_minor_size)) +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.8, hjust = 0.8)) +
  theme(legend.title = element_text(size = text_major_size)) +
  theme(legend.text = element_text(size = text_minor_size)) +
  theme(strip.text = element_text(size = text_minor_size, color = '#000000')) +
  theme(strip.background = element_rect(fill = '#dddddd')) +
  facet_grid(rows = vars(length_factor)) +
  xlab('Organism size') +
  ylab('Number of finished replicates')

# Plot the evolved restraint buffer as boxplots
  # x-axis = organism size
  # y-axis = average evolved restraint buffer for each replicate
  # facet rows = genome length
ggplot(df2[df2$generation == 10000,], aes(x = size_factor, y = restraint_value)) +
  geom_hline(aes(yintercept = 0), alpha = 0.5, linetype = 'dashed') +
  geom_boxplot(aes(fill = size_factor)) +
  xlab('Organism size') +
  ylab('Evolved restraint buffer') +
  scale_fill_manual(values = color_map) +
  labs(fill = 'Organism size') +
  theme_light() +
  theme(legend.position = 'none') +
  theme(panel.grid.major.x = element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(axis.title = element_text(size = text_major_size)) +
  theme(axis.text = element_text(size = text_minor_size)) +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.8, hjust = 0.8)) +
  theme(legend.title = element_text(size = text_major_size)) +
  theme(legend.text = element_text(size = text_minor_size)) +
  theme(strip.text = element_text(size = text_minor_size, color = '#000000')) +
  theme(strip.background = element_rect(fill = '#dddddd')) +
  facet_grid(vars(length_factor)) 
ggsave('./evolution/plots/genome_length_boxplot_facet_length.png', units = 'in', width = 6, height = 10) 
ggsave('./evolution/plots/genome_length_boxplot_facet_length.pdf', units = 'in', width = 6, height = 10)
# Same plot but y axes are allowed to vary between rows
ggplot(df2[df2$generation == 10000,], aes(x = size_factor, y = restraint_value)) +
  geom_hline(aes(yintercept = 0), alpha = 0.5, linetype = 'dashed') +
  geom_boxplot(aes(fill = size_factor)) +
  xlab('Organism size') +
  ylab('Evolved restraint buffer') +
  scale_fill_manual(values = color_map) +
  labs(fill = 'Organism size') +
  theme_light() +
  theme(legend.position = 'none') +
  theme(panel.grid.major.x = element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(axis.title = element_text(size = text_major_size)) +
  theme(axis.text = element_text(size = text_minor_size)) +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.8, hjust = 0.8)) +
  theme(legend.title = element_text(size = text_major_size)) +
  theme(legend.text = element_text(size = text_minor_size)) +
  theme(strip.text = element_text(size = text_minor_size, color = '#000000')) +
  theme(strip.background = element_rect(fill = '#dddddd')) +
  facet_grid(vars(length_factor), scales = 'free_y') 
ggsave('./evolution/plots/genome_length_boxplot_facet_length_free_y.png', units = 'in', width = 6, height = 10) 
ggsave('./evolution/plots/genome_length_boxplot_facet_length_free_y.pdf', units = 'in', width = 6, height = 10)

# Plot the evolved restraint buffer as boxplots
  # x-axis = genome length
  # y-axis = average evolved restraint buffer for each replicate
  # facet rows = organism size
ggplot(df2[df2$generation == 10000,], aes(x = length_factor, y = restraint_value)) +
  geom_hline(aes(yintercept = 0), alpha = 0.5, linetype = 'dashed') +
  geom_boxplot(aes(fill = size_factor)) +
  xlab('Genome length') +
  ylab('Evolved restraint buffer') +
  scale_fill_manual(values = color_map) +
  labs(fill = 'Organism size') +
  theme_light() +
  theme(legend.position = 'none') +
  theme(panel.grid.major.x = element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(axis.title = element_text(size = text_major_size)) +
  theme(axis.text = element_text(size = text_minor_size)) +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.8, hjust = 0.8)) +
  theme(legend.title = element_text(size = text_major_size)) +
  theme(legend.text = element_text(size = text_minor_size)) +
  theme(strip.text = element_text(size = text_minor_size, color = '#000000')) +
  theme(strip.background = element_rect(fill = '#dddddd')) +
  facet_grid(vars(size_factor)) 
ggsave('./evolution/plots/genome_length_boxplot_facet_size.png', units = 'in', width = 6, height = 10) 
ggsave('./evolution/plots/genome_length_boxplot_facet_size.pdf', units = 'in', width = 6, height = 10)
# Same plot, but allow each row's y-axis to vary
ggplot(df2[df2$generation == 10000,], aes(x = length_factor, y = restraint_value)) +
  geom_hline(aes(yintercept = 0), alpha = 0.5, linetype = 'dashed') +
  geom_boxplot(aes(fill = size_factor)) +
  xlab('Genome length') +
  ylab('Evolved restraint buffer') +
  scale_fill_manual(values = color_map) +
  labs(fill = 'Organism size') +
  theme_light() +
  theme(legend.position = 'none') +
  theme(panel.grid.major.x = element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(axis.title = element_text(size = text_major_size)) +
  theme(axis.text = element_text(size = text_minor_size)) +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.8, hjust = 0.8)) +
  theme(legend.title = element_text(size = text_major_size)) +
  theme(legend.text = element_text(size = text_minor_size)) +
  theme(strip.text = element_text(size = text_minor_size, color = '#000000')) +
  theme(strip.background = element_rect(fill = '#dddddd')) +
  facet_grid(vars(size_factor), scales = 'free_y') 
ggsave('./evolution/plots/genome_length_boxplot_facet_size_free_y.png', units = 'in', width = 6, height = 10) 
ggsave('./evolution/plots/genome_length_boxplot_facet_size_free_y.pdf', units = 'in', width = 6, height = 10)

# Plot the evolved restraint buffer as boxplots, one file per genome length
  # x-axis = organism size
  # y-axis = average evolved restraint buffer for each replicate
for(genome_length in unique(df2$LENGTH)){
  ggplot(df2[df2$generation == 10000 & df2$LENGTH == genome_length,], aes(x = size_factor, y = restraint_value)) +
    geom_hline(aes(yintercept = 0), alpha = 0.5, linetype = 'dashed') +
    geom_boxplot(aes(fill = size_factor)) +
    xlab('Organism size') +
    ylab('Evolved restraint buffer') +
    scale_fill_manual(values = color_map) +
    labs(fill = 'Organism size') +
    theme_light() +
    theme(legend.position = 'none') +
    theme(panel.grid.major.x = element_blank()) +
    theme(panel.grid.minor.x = element_blank()) +
    theme(axis.title = element_text(size = text_major_size)) +
    theme(axis.text = element_text(size = text_minor_size)) +
    theme(axis.text.x = element_text(angle = 60, vjust = 0.8, hjust = 0.8)) +
    theme(legend.title = element_text(size = text_major_size)) +
    theme(legend.text = element_text(size = text_minor_size)) +
    theme(strip.text = element_text(size = text_minor_size, color = '#000000')) +
    theme(strip.background = element_rect(fill = '#dddddd')) +
    facet_grid(vars(length_factor)) 
  ggsave(paste0('./evolution/plots/genome_length_boxplot_single_length_', genome_length, '.png'), units = 'in', width = 6, height = 6) 
  ggsave(paste0('./evolution/plots/genome_length_boxplot_single_length_', genome_length, '.pdf'), units = 'in', width = 6, height = 6) 
}

# Plot the evolved restraint buffer as boxplots, one file per organism size
  # x-axis = genome length
  # y-axis = average evolved restraint buffer for each replicate
for(org_size in unique(df2$MCSIZE)){
  ggplot(df2[df2$generation == 10000 & df2$MCSIZE == org_size,], aes(x = length_factor, y = restraint_value)) +
    geom_hline(aes(yintercept = 0), alpha = 0.5, linetype = 'dashed') +
    geom_boxplot(aes(fill = size_factor)) +
    xlab('Genome length') +
    ylab('Evolved restraint buffer') +
    scale_fill_manual(values = color_map) +
    labs(fill = 'Organism size') +
    theme_light() +
    theme(legend.position = 'none') +
    theme(panel.grid.major.x = element_blank()) +
    theme(panel.grid.minor.x = element_blank()) +
    theme(axis.title = element_text(size = text_major_size)) +
    theme(axis.text = element_text(size = text_minor_size)) +
    theme(axis.text.x = element_text(angle = 60, vjust = 0.8, hjust = 0.8)) +
    theme(legend.title = element_text(size = text_major_size)) +
    theme(legend.text = element_text(size = text_minor_size)) +
    theme(strip.text = element_text(size = text_minor_size, color = '#000000')) +
    theme(strip.background = element_rect(fill = '#dddddd')) +
    facet_grid(vars(size_factor)) 
  ggsave(paste0('./evolution/plots/genome_length_boxplot_single_size_', org_size, '.png'), units = 'in', width = 6, height = 6) 
  ggsave(paste0('./evolution/plots/genome_length_boxplot_single_size_', org_size, '.pdf'), units = 'in', width = 6, height = 6)
}
