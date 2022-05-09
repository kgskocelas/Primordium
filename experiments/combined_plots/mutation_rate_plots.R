rm(list = ls())

# Load relevant libraries
library(dplyr)
library(ggplot2)
library(ggridges)
library(scales)
library(khroma)
library(cowplot)

####################################
############# GERM #################
####################################

# Load the data
df_germ = read.csv('../2021_02_16__germ_mut_fin/evolution/data/scraped_evolution_data.csv')

# Create directory to save plots
if(!dir.exists('./plots')){
  dir.create('./plots')
}

# Trim off NAs (artifiacts of how we scraped the data) and trim to only have gen 10,000
df_germ2 = df_germ[!is.na(df_germ$MCSIZE) & df_germ$generation == 10000,]
# Ignore data for size 8x8 and 1024x1024
df_germ2 = df_germ2[df_germ2$MCSIZE != 8 & df_germ2$MCSIZE != 1024,]

## Set variables to make plotting easier
# Calculate restraint value (x - 60 because genome length is 100 here)
df_germ2$restraint_value = df_germ2$ave_ones - 60
# Make a nice, clean factor for size
df_germ2$size_str = paste0(df_germ2$MCSIZE, 'x', df_germ2$MCSIZE)
df_germ2$size_factor = factor(df_germ2$size_str, levels = c('16x16', '32x32', '64x64', '128x128', '256x256', '512x512', '1024x1024'))
df_germ2$size_factor_reversed = factor(df_germ2$size_str, levels = rev(c('16x16', '32x32', '64x64', '128x128', '256x256', '512x512', '1024x1024')))
df_germ2$germ_mut_str = paste('GERM MUT', df_germ2$MUT)
df_germ2$mut_factor = factor(df_germ2$MUT, levels = c(0.01, 0.02, 0.05, 0.10, 0.20, 0.50, 1.00))
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
boxplot_color = '#9ecae1'

ggp_germ = ggplot(df_germ2[df_germ2$generation == 10000 & df_germ2$MCSIZE == 256,], aes(x = mut_factor, y = restraint_value)) +
  geom_hline(aes(yintercept = 0), alpha = 0.5, linetype = 'dashed') +
  geom_boxplot(aes(fill = size_factor)) +
  xlab('Germ mutation rate') +
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
  scale_y_continuous(limits = c(-22, 22))

####################################
############# SOMA #################
####################################

# Load the data
df_soma = read.csv('../2021_02_27__soma_mut_fin/evolution/data/scraped_evolution_data.csv')

# Trim off NAs (artifiacts of how we scraped the data) and trim to only have gen 10,000
df_soma2 = df_soma[!is.na(df_soma$MCSIZE) & df_soma$generation == 10000,]
# Ignore data for size 8x8 and 1024x1024
df_soma2 = df_soma2[df_soma2$MCSIZE != 8 & df_soma2$MCSIZE != 1024,]


## Set variables to make plotting easier
# Calculate restraint value (x - 60 because genome length is 100 here)
df_soma2$restraint_value = df_soma2$ave_ones - 60
# Make a nice, clean factor for size
df_soma2$size_str = paste0(df_soma2$MCSIZE, 'x', df_soma2$MCSIZE)
df_soma2$size_factor = factor(df_soma2$size_str, levels = c('16x16', '32x32', '64x64', '128x128', '256x256', '512x512', '1024x1024'))
df_soma2$size_factor_reversed = factor(df_soma2$size_str, levels = rev(c('16x16', '32x32', '64x64', '128x128', '256x256', '512x512', '1024x1024')))
df_soma2$soma_mut_str = paste('soma CELLMUT', df_soma2$CELLMUT)
df_soma2$mut_factor = factor(df_soma2$CELLMUT, levels = c(0.01, 0.02, 0.05, 0.10, 0.20, 0.50, 1.00))

# x-axis = soma mutation rate
# y-axis = average evolved restraint buffer for each replicate
ggp_soma = ggplot(df_soma2[df_soma2$generation == 10000 & df_soma2$MCSIZE == 256,], aes(x = mut_factor, y = restraint_value)) +
  geom_hline(aes(yintercept = 0), alpha = 0.5, linetype = 'dashed') +
  geom_boxplot(aes(fill = size_factor)) +
  xlab('Soma mutation rate') +
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
  scale_y_continuous(limits = c(-22, 22))


####################################
########## COMBINED ################
####################################

plot_grid(ggp_germ, ggp_soma + theme(axis.title.y = element_blank()), labels = c('A', 'B'), label_size = text_major_size, nrow = 1)
ggsave('./plots/combined_mut_rate_plots.png', units = 'in', width = 6, height = 6) 
ggsave('./plots/combined_mut_rate_plots.pdf', units = 'in', width = 6, height = 6) 

