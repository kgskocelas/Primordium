rm(list = ls())

# Load libraries
library(ggplot2)
library(cowplot)
library(khroma)
library(ggsignif)


size_vec = c(16, 32, 64, 128, 256, 512)
num_samples = 3000
num_batches = 30 
batch_size = num_samples / num_batches
y_position_even = 80 
y_position_odd = 85
max_gen = 10000

########################################
############ FINITE POP ################
########################################

# Load data
df_finite = NA
for(size in size_vec){
  for(batch_id in 1:num_batches){
    if(!is.data.frame(df_finite)){
      df_finite = read.csv(paste0('../2021_03_07__inf_population/simulation/data/size_', size, '__batch_', batch_id, '.csv'))
      df_finite$batch_id = batch_id
    } else{
      df_finite_tmp = read.csv(paste0('../2021_03_07__inf_population/simulation/data/size_', size, '__batch_', batch_id, '.csv'))
      df_finite_tmp$batch_id = batch_id
      df_finite = rbind(df_finite, df_finite_tmp)
    }
  }
}

# Summarize data, grabbing the mean from the simulation
df_finite_summary = data.frame(data = matrix(nrow = 0, ncol = 4))
colnames(df_finite_summary) = c('mc_size', 'generation', 'mean_ones', 'batch_id')
for(size in size_vec){
  size_mask = df_finite$mc_size == size
  for(batch_id in unique(df_finite$batch_id)){
    mask = size_mask & df_finite$batch_id == batch_id 
    cat(paste0('(', batch_id, ') Size: ', size, '\n\t'))
    for(gen in unique(df_finite$generation)){
      df_finite_summary[nrow(df_finite_summary) + 1,] = c(size, gen, weighted.mean(0:100, df_finite[mask & df_finite$generation == gen,]$frac_of_pop), batch_id)
      cat(gen, ' ')
    }
    cat('\n')
  }
}

# Trim to final generation
df_finite_summary = df_finite_summary[df_finite_summary$generation == max_gen,]

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

# Calculate stats
df_finite_wilcox = data.frame(data = matrix(nrow = 0, ncol = 7))
adjacent_counter = 0
colnames(df_finite_wilcox) = c('size_a', 'size_b', 'p_value_corrected', 'p_value_raw', 'W', 'is_adjacent', 'y_position')
for(size_idx_a in 1:(length(size_vec) - 1)){
  size_a = size_vec[size_idx_a]
  size_a_str = paste0(size_a, 'x', size_a)
  for(size_idx_b in (size_idx_a + 1):length(size_vec)){
    size_b = size_vec[size_idx_b]
    size_b_str = paste0(size_b, 'x', size_b)
    res = wilcox.test(
      df_finite_summary[df_finite_summary$mc_size == size_a,]$restraint_value, 
      df_finite_summary[df_finite_summary$mc_size == size_b,]$restraint_value, 
      alternative = 'two.sided') 
    is_adjacent = (size_b == size_a * 2)
    y_position = y_position_even
    if(is_adjacent == T){
      if(adjacent_counter %% 2 != 0){
        y_position = y_position_odd
      }
      adjacent_counter = adjacent_counter + 1
    }
    df_finite_wilcox[nrow(df_finite_wilcox) + 1,] = c(size_a_str, size_b_str, 0, res$p.value, as.numeric(res$statistic)[1], is_adjacent, as.numeric(y_position))
  }
}
df_finite_wilcox$p_value_corrected = p.adjust(df_finite_wilcox$p_value_raw, method = 'holm')
df_finite_wilcox$less_0.01 = df_finite_wilcox$p_value_corrected < 0.01
df_finite_wilcox$label = 'ns'
df_finite_wilcox[df_finite_wilcox$p_value_corrected <= 0.05,]$label = '*'
df_finite_wilcox[df_finite_wilcox$p_value_corrected <= 0.01,]$label = '**'
df_finite_wilcox[df_finite_wilcox$p_value_corrected <= 0.001,]$label = '***'
df_finite_wilcox$size_a_factor = factor(df_finite_wilcox$size_a, levels = c('16x16', '32x32', '64x64', '128x128', '256x256', '512x512'))
df_finite_wilcox$size_b_factor = factor(df_finite_wilcox$size_b, levels = c('16x16', '32x32', '64x64', '128x128', '256x256', '512x512'))
print(df_finite_wilcox)

ggp_finite = ggplot(df_finite_summary[df_finite_summary$generation == 10000,], aes(x = size_factor, y = restraint_value ))  +
  geom_hline(aes(yintercept = 0), alpha = 0.5, linetype = 'dashed') + 
  geom_boxplot(aes(fill = size_factor)) + 
  #geom_text(aes(x = size_factor, y = restraint_value + 2, label = round(restraint_value, 1))) + 
  geom_signif(
    data = df_finite_wilcox[df_finite_wilcox$is_adjacent == T,],
    aes(annotations = label, xmin=size_a_factor, xmax=size_b_factor,y_position=as.numeric(y_position)),
    manual=T,
    inherit.aes=FALSE
  ) +
  xlab('Organism size') + 
  ylab('Average restraint buffer') + 
  labs(color = 'Organism size') + 
  scale_fill_manual(values = color_map) + 
  scale_y_continuous(limits = c(0, 85)) +
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

ggp_finite

########################################
########### INFINITE POP ###############
########################################

# Load data
df_infinite = NA
for(size in size_vec){
  for(batch_id in 1:num_batches){
    if(!is.data.frame(df_infinite)){
      df_infinite = read.csv(paste0('../2021_03_08__inf_genome_inf_pop/simulation/data/size_', size, '__batch_', batch_id, '.csv'))
      df_infinite$batch_id = batch_id
    } else{
      df_infinite_tmp = read.csv(paste0('../2021_03_08__inf_genome_inf_pop/simulation/data/size_', size, '__batch_', batch_id, '.csv'))
      df_infinite_tmp$batch_id = batch_id
      df_infinite = rbind(df_infinite, df_infinite_tmp)
    }
  }
}

# Summarize data, grabbing the mean from the simulation
df_infinite_summary = data.frame(data = matrix(nrow = 0, ncol = 4))
colnames(df_infinite_summary) = c('mc_size', 'generation', 'mean_ones', 'batch_id')
for(size in size_vec){
  size_mask = df_infinite$mc_size == size
  for(batch_id in unique(df_infinite$batch_id)){
    mask = size_mask & df_infinite$batch_id == batch_id
    cat(paste0('(', batch_id, ') Size: ', size, '\n\t'))
    for(gen in unique(df_infinite$generation)){
      df_infinite_summary[nrow(df_infinite_summary) + 1,] = c(size, gen, weighted.mean(-100:550, df_infinite[mask & df_infinite$generation == gen,]$frac_of_pop), batch_id)
      cat(gen, ' ')
    }
    cat('\n')
  }
}

# Trim to final generation
df_infinite_summary = df_infinite_summary[df_infinite_summary$generation == max_gen,]

# Create variables to make plotting easier
df_infinite_summary$restraint_value = df_infinite_summary$mean_ones 
df_infinite_summary$size_str = paste0(df_infinite_summary$mc_size, 'x', df_infinite_summary$mc_size)
df_infinite_summary$size_factor = factor(df_infinite_summary$size_str, levels = paste0(size_vec, 'x', size_vec))

# Calculate stats
df_infinite_wilcox = data.frame(data = matrix(nrow = 0, ncol = 7))
adjacent_counter = 0
colnames(df_infinite_wilcox) = c('size_a', 'size_b', 'p_value_corrected', 'p_value_raw', 'W', 'is_adjacent', 'y_position')
for(size_idx_a in 1:(length(size_vec) - 1)){
  size_a = size_vec[size_idx_a]
  size_a_str = paste0(size_a, 'x', size_a)
  for(size_idx_b in (size_idx_a + 1):length(size_vec)){
    size_b = size_vec[size_idx_b]
    size_b_str = paste0(size_b, 'x', size_b)
    res = wilcox.test(
      df_infinite_summary[df_infinite_summary$mc_size == size_a,]$restraint_value, 
      df_infinite_summary[df_infinite_summary$mc_size == size_b,]$restraint_value, 
      alternative = 'two.sided') 
    is_adjacent = (size_b == size_a * 2)
    y_position = y_position_even
    if(is_adjacent == T){
      if(adjacent_counter %% 2 != 0){
        y_position = y_position_odd
      }
      adjacent_counter = adjacent_counter + 1
    }
    df_infinite_wilcox[nrow(df_infinite_wilcox) + 1,] = c(size_a_str, size_b_str, 0, res$p.value, as.numeric(res$statistic)[1], is_adjacent, as.numeric(y_position))
  }
}
df_infinite_wilcox$p_value_corrected = p.adjust(df_infinite_wilcox$p_value_raw, method = 'holm')
df_infinite_wilcox$less_0.01 = df_infinite_wilcox$p_value_corrected < 0.01
df_infinite_wilcox$label = 'ns'
df_infinite_wilcox[df_infinite_wilcox$p_value_corrected <= 0.05,]$label = '*'
df_infinite_wilcox[df_infinite_wilcox$p_value_corrected <= 0.01,]$label = '**'
df_infinite_wilcox[df_infinite_wilcox$p_value_corrected <= 0.001,]$label = '***'
df_infinite_wilcox$size_a_factor = factor(df_infinite_wilcox$size_a, levels = c('16x16', '32x32', '64x64', '128x128', '256x256', '512x512'))
df_infinite_wilcox$size_b_factor = factor(df_infinite_wilcox$size_b, levels = c('16x16', '32x32', '64x64', '128x128', '256x256', '512x512'))
print(df_infinite_wilcox)

ggp_infinite = ggplot(df_infinite_summary[df_infinite_summary$generation == 10000,], aes(x = size_factor, y = restraint_value ))  +
  geom_hline(aes(yintercept = 0), alpha = 0.5, linetype = 'dashed') + 
  geom_boxplot(aes(fill = size_factor)) + 
  #geom_text(aes(x = size_factor, y = restraint_value + 2, label = round(restraint_value, 1))) + 
  geom_signif(
    data = df_infinite_wilcox[df_infinite_wilcox$is_adjacent == T,],
    aes(annotations = label, xmin=size_a_factor, xmax=size_b_factor,y_position=as.numeric(y_position)),
    manual=T,
    inherit.aes=FALSE
  ) +
  xlab('Organism size') + 
  ylab('Average restraint buffer') + 
  labs(color = 'Organism size') + 
  scale_fill_manual(values = color_map) + 
  scale_y_continuous(limits = c(0, 85)) +
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
ggp_infinite

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
