rm(list = ls())

# Load relevant libraries
library(dplyr)
library(ggplot2)
library(ggridges)
library(scales)
library(khroma)
library(ggsignif)

# Load the data
df = read.csv('evolution/data/scraped_evolution_data.csv')

# Create directory to save plots
if(!dir.exists('./evolution/plots')){
  dir.create('./evolution/plots')
}

# Trim off NAs (artifiacts of how we scraped the data) and trim to only have gen 10,000
df2 = df[!is.na(df$MCSIZE) & df$generation == 10000,]
# Ignore data for size 8x8 and 1024x1024
df2 = df2[df2$MCSIZE != 8 & df2$MCSIZE != 1024,]

# Group the data by size and summarize
data_grouped = dplyr::group_by(df2, MCSIZE)
data_summary = dplyr::summarize(data_grouped, mean_ones = mean(ave_ones), n = dplyr::n())

# Print any data that is missing!
if(sum(data_summary$n != 100) != 0){
  cat('Error! Missing data detected!\n')
  print(data_summary[data_summary$n != 100,])
}

## Set variables to make plotting easier
# Calculate restraint value (x - 60 because genome length is 100 here)
df2$restraint_value = df2$ave_ones - 60
# Make a nice, clean factor for size
df2$size_str = paste0(df2$MCSIZE, 'x', df2$MCSIZE)
df2$size_factor = factor(df2$size_str, levels = c('16x16', '32x32', '64x64', '128x128', '256x256', '512x512', '1024x1024'))
df2$size_factor_reversed = factor(df2$size_str, levels = rev(c('16x16', '32x32', '64x64', '128x128', '256x256', '512x512', '1024x1024')))
data_summary$size_str = paste0(data_summary$MCSIZE, 'x', data_summary$MCSIZE)
data_summary$size_factor = factor(data_summary$size_str, levels = c('16x16', '32x32', '64x64', '128x128', '256x256', '512x512', '1024x1024'))
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
y_position_even = 24
y_position_odd = 22
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

# Plot the number of replicates for each organism size
ggplot(data_summary, aes(x = size_factor, y = n)) +
  geom_col(aes(fill = size_factor)) +
  geom_text(aes(y = n + 2, label = n)) +
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
  xlab('Organism size') +
  ylab('Number of finished replicates')


# Plot the evolved restraint buffer for all reps of each org size as boxplots
ggplot(df2[df2$generation == 10000,], aes(x = size_factor, y = restraint_value)) +
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
  theme(strip.background = element_rect(fill = '#dddddd'))  
  #ggtitle('Finite Genome Avg Restraint by MC Size') 
ggsave('./evolution/plots/finite_extended_size_boxplot.png', units = 'in', width = 6, height = 6) 
ggsave('./evolution/plots/finite_extended_size_boxplot.pdf', units = 'in', width = 6, height = 6)

# Plot same data as a raincloud plot
ggplot(df2[df2$generation == 10000,], aes(x = restraint_value, y = size_factor_reversed)) +
  geom_vline(aes(xintercept = 0), alpha = 0.5, linetype = 'dashed') +
  geom_density_ridges2(aes(fill = size_factor), scale = 0.5, position = position_nudge(y = 0.3)) +
  geom_jitter(aes(color = size_factor), height = 0.18, size = 0.3) + 
  xlab('Evolved restraint buffer') +
  ylab('Organism size') +
  scale_fill_manual(values = color_map) +
  scale_color_manual(values = color_map) +
  labs(fill = 'Organism size') +
  theme_light() +
  theme(legend.position = 'none') +
  theme(panel.grid.major.x = element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(axis.title = element_text(size = text_major_size)) +
  theme(axis.text = element_text(size = text_minor_size)) +
  theme(legend.title = element_text(size = text_major_size)) +
  theme(legend.text = element_text(size = text_minor_size)) +
  theme(strip.text = element_text(size = text_minor_size, color = '#000000')) +
  theme(strip.background = element_rect(fill = '#dddddd'))  
ggsave('./evolution/plots/finite_extended_size_raincloud.png', units = 'in', width = 6, height = 6) 
ggsave('./evolution/plots/finite_extended_size_raincloud.pdf', units = 'in', width = 6, height = 6)

