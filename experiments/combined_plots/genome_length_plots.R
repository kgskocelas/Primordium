rm(list = ls())

# Load relevant libraries
library(dplyr)
library(ggplot2)
library(ggridges)
library(scales)
library(khroma)
library(cowplot)
library(ggsignif)

# Load the data
df = read.csv('../2021_02_27__genome_length/evolution/data/scraped_evolution_data_length_50.csv')
df = rbind(df, read.csv('../2021_02_27__genome_length/evolution/data/scraped_evolution_data_length_200.csv'))
df = rbind(df, read.csv('../2021_02_27__genome_length/evolution/data/scraped_evolution_data_length_100.csv'))
df = rbind(df, read.csv('../2021_02_27__genome_length/evolution/data/scraped_evolution_data_length_400.csv'))
df = rbind(df, read.csv('../2021_02_27__genome_length/evolution/data/scraped_evolution_data_length_25.csv'))

# Create directory to save plots
if(!dir.exists('./plots')){
  dir.create('./plots')
}

# Trim off NAs (artifiacts of how we scraped the data) and trim to only have gen 10,000
df2 = df[!is.na(df$MCSIZE) & df$generation == 10000,]
# Ignore data for size 8x8 and 1024x1024
df2 = df2[df2$MCSIZE != 8 & df2$MCSIZE != 1024,]

## Set variables to make plotting easier
# Calculate restraint value (x - 60% of the genome length)
df2$restraint_value = df2$ave_ones - (df2$LENGTH * 0.6)
# Make a nice, clean factor for size
df2$size_str = paste0(df2$MCSIZE, 'x', df2$MCSIZE)
df2$size_factor = factor(df2$size_str, levels = c('16x16', '32x32', '64x64', '128x128', '256x256', '512x512', '1024x1024'))
df2$size_factor_reversed = factor(df2$size_str, levels = rev(c('16x16', '32x32', '64x64', '128x128', '256x256', '512x512', '1024x1024')))
df2$length_str = paste0(df2$LENGTH, '-bit')
df2$length_factor = factor(df2$length_str, levels = c('25-bit', '50-bit', '100-bit', '200-bit', '400-bit'))
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

# Run the stats necessary for these plots (rest are in supplement)
size_vec = c(16, 32, 64, 128, 256, 512)
y_position_even = 35
y_position_odd = 40
adjacent_counter = 0
df_test = df2
df_wilcox = data.frame(data = matrix(nrow = 0, ncol = 7))
colnames(df_wilcox) = c('size_a', 'size_b', 'p_value_corrected', 'p_value_raw', 'W', 'is_adjacent', 'y_position')
for(size_idx_a in 1:(length(size_vec) - 1)){
  size_a = size_vec[size_idx_a]
  size_a_str = paste0(size_a, 'x', size_a)
  for(size_idx_b in (size_idx_a + 1):length(size_vec)){
    size_b = size_vec[size_idx_b]
    size_b_str = paste0(size_b, 'x', size_b)
    res = wilcox.test(df_test[df_test$MCSIZE == size_a,]$restraint_value, df_test[df_test$MCSIZE == size_b,]$restraint_value, alternative = 'two.sided') 
    is_adjacent = (size_b == size_a * 2)
    y_position = y_position_even
    if(is_adjacent == T){
      if(adjacent_counter %% 2 != 0){
        y_position = y_position_odd
      }
      adjacent_counter = adjacent_counter + 1
    }
    df_wilcox[nrow(df_wilcox) + 1,] = c(size_a_str, size_b_str, 0, res$p.value, as.numeric(res$statistic)[1], is_adjacent, as.numeric(y_position))
  }
}
df_wilcox$p_value_corrected = p.adjust(df_wilcox$p_value_raw, method = 'holm')
df_wilcox$less_0.01 = df_wilcox$p_value_corrected < 0.01
df_wilcox$label = 'ns'
df_wilcox[df_wilcox$p_value_corrected <= 0.05,]$label = '*'
df_wilcox[df_wilcox$p_value_corrected <= 0.01,]$label = '**'
df_wilcox[df_wilcox$p_value_corrected <= 0.001,]$label = '***'
df_wilcox$size_a_factor = factor(df_wilcox$size_a, levels = c('16x16', '32x32', '64x64', '128x128', '256x256', '512x512'))
df_wilcox$size_b_factor = factor(df_wilcox$size_b, levels = c('16x16', '32x32', '64x64', '128x128', '256x256', '512x512'))
print(df_wilcox)

# x-axis = organism size
# y-axis = average evolved restraint buffer for each replicate
ggp_400_bit = ggplot(df2[df2$generation == 10000 & df2$LENGTH == 400,], aes(x = size_factor, y = restraint_value)) +
  geom_hline(aes(yintercept = 0), alpha = 0.5, linetype = 'dashed') +
  geom_boxplot(aes(fill = size_factor)) +
  geom_signif(
    data = df_wilcox[df_wilcox$is_adjacent == T,],
    aes(annotations = label, xmin=size_a_factor, xmax=size_b_factor,y_position=as.numeric(y_position)),
    manual=T,
    inherit.aes=FALSE
  ) +
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
  scale_y_continuous(limits=c(-42, 40))
ggp_400_bit


# Run the stats necessary for these plots (rest are in supplement)
size_vec = c(16, 32, 64, 128, 256, 512)
y_position_even = 35
y_position_odd = 40
adjacent_counter = 0
df_test = df2
df_wilcox = data.frame(data = matrix(nrow = 0, ncol = 7))
colnames(df_wilcox) = c('size_a', 'size_b', 'p_value_corrected', 'p_value_raw', 'W', 'is_adjacent', 'y_position')
for(size_idx_a in 1:(length(size_vec) - 1)){
  size_a = size_vec[size_idx_a]
  size_a_str = paste0(size_a, 'x', size_a)
  for(size_idx_b in (size_idx_a + 1):length(size_vec)){
    size_b = size_vec[size_idx_b]
    size_b_str = paste0(size_b, 'x', size_b)
    res = wilcox.test(df_test[df_test$MCSIZE == size_a,]$restraint_value, df_test[df_test$MCSIZE == size_b,]$restraint_value, alternative = 'two.sided') 
    is_adjacent = (size_b == size_a * 2)
    y_position = y_position_even
    if(is_adjacent == T){
      if(adjacent_counter %% 2 != 0){
        y_position = y_position_odd
      }
      adjacent_counter = adjacent_counter + 1
    }
    df_wilcox[nrow(df_wilcox) + 1,] = c(size_a_str, size_b_str, 0, res$p.value, as.numeric(res$statistic)[1], is_adjacent, as.numeric(y_position))
  }
}
df_wilcox$p_value_corrected = p.adjust(df_wilcox$p_value_raw, method = 'holm')
df_wilcox$less_0.01 = df_wilcox$p_value_corrected < 0.01
df_wilcox$label = 'ns'
df_wilcox[df_wilcox$p_value_corrected <= 0.05,]$label = '*'
df_wilcox[df_wilcox$p_value_corrected <= 0.01,]$label = '**'
df_wilcox[df_wilcox$p_value_corrected <= 0.001,]$label = '***'
df_wilcox$size_a_factor = factor(df_wilcox$size_a, levels = c('16x16', '32x32', '64x64', '128x128', '256x256', '512x512'))
df_wilcox$size_b_factor = factor(df_wilcox$size_b, levels = c('16x16', '32x32', '64x64', '128x128', '256x256', '512x512'))
print(df_wilcox)

# x-axis = genome length
# y-axis = average evolved restraint buffer for each replicate
ggp_size_256 = ggplot(df2[df2$generation == 10000 & df2$MCSIZE == 256,], aes(x = length_factor, y = restraint_value)) +
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
  scale_y_continuous(limits=c(-42, 40))
ggp_size_256

####################################
########## COMBINED ################
####################################

plot_grid(ggp_size_256, ggp_400_bit + theme(axis.title.y = element_blank()), labels = c('A', 'B'), label_size = text_major_size, nrow = 1)
ggsave('./plots/combined_genome_length_plots.png', units = 'in', width = 6, height = 6) 
ggsave('./plots/combined_genome_length_plots.pdf', units = 'in', width = 6, height = 6) 
