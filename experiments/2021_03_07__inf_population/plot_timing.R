rm(list = ls())
library(ggplot2)
library(dplyr)

# Load in data
df_summary = read.csv('./timing/data/timing_data_summary.csv')

# Create a ffew variables to make plotting easier
df_summary$size_str = paste0(df_summary$cells_side, 'x', df_summary$cells_side)
df_summary$size_factor = factor(df_summary$size_str, levels = c('8x8', '16x16', '32x32', '64x64', '128x128', '256x256', '512x512'))
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
text_major_size = 18
text_minor_size = 16 

# Create directory to save plots
if(!dir.exists('./timing')){
  dir.create('./timing')
}
if(!dir.exists('./timing/plots')){
  dir.create('./timing/plots')
}

ggplot(df_summary, aes(x = ancestor_1s - 60, y = mean_rep_rate, color = as.factor(size_factor))) + 
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

ggplot(df_summary[df_summary$ancestor_1s >= 40,], aes(x = ancestor_1s - 60, y = rep_rate_diff, color = as.factor(size_factor))) +
  geom_vline(aes(xintercept = 0), alpha = 0.5, linetype = 'dashed') +
  geom_line(size = 1.0) + 
  xlab('Simulated restraint bufffer') + 
  ylab('Change in fitness') +
  labs(color = 'Organism\n  size') +
  scale_color_manual(values = color_map) +
  theme_light() +
  theme(panel.grid.major.x = element_blank()) +
  theme(axis.title = element_text(size = text_major_size)) +
  theme(axis.text = element_text(size = text_minor_size)) +
  theme(legend.title = element_text(size = text_major_size)) +
  theme(legend.text = element_text(size = text_minor_size)) +
  theme(legend.position = 'bottom') +
  guides(color = guide_legend(nrow = 3)) 
ggsave('timing/plots/replication_rate_diff_trimmed.png', units = 'in', width = 6, height = 6) 
ggsave('timing/plots/replication_rate_diff_trimmed.pdf', units = 'in', width = 6, height = 6)

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
  #theme(legend.position = 'bottom') +
  theme(legend.position = 'none') 
  #guides(color = guide_legend(nrow = 3)) 
ggsave('timing/plots/replication_rate_diff_trimmed_two_col.png', units = 'in', width = 6, height = 6) 
ggsave('timing/plots/replication_rate_diff_trimmed_two_col.pdf', units = 'in', width = 6, height = 6)
