rm(list = ls())

library(ggplot2)

#rgb(r, g, b, maxColorValue = 255)
threshold = 60

df = data.frame(data = matrix(nrow = 101, ncol = 0))
df$num_ones = 0:100
df$red = 0
df$green = 0
df$blue = 0
df$color = rgb(0,0,0)
for(ones in 0:100){
  cat(ones, ' ')
  if(ones < threshold){
    df[ones + 1,]$red = 255 - (threshold - ones) * (255 / threshold)
    df[ones + 1,]$green = 0
    df[ones + 1,]$blue = 0 + (threshold - ones) *   (100 / threshold)
  }
  else{
    df[ones + 1,]$red = 255 - (ones - threshold) * (255 / (100 - threshold))
    df[ones + 1,]$green = 255 - (ones - threshold) * (255 / (100 - threshold))
    df[ones + 1,]$blue = 255 - (ones - threshold) * (255 / (100 - threshold))
  }
  df[ones + 1,]$color = rgb(
    df[ones + 1,]$red,
    df[ones + 1,]$green,
    df[ones + 1,]$blue, 
    maxColorValue = 255
  )
}

# Full range based on number of ones
ggplot(df, aes(x = num_ones, y = 1, fill = as.factor(num_ones))) + 
  geom_tile() + 
  scale_fill_manual(values = df$color) +
  scale_x_continuous(expand = c(0,0), breaks = seq(0,100,10)) +
  scale_y_continuous(expand = c(0,0)) +
  xlab('Number of ones') +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(legend.position = 'none') 
ggsave('./misc_plots/num_ones_colormap_full.pdf', units = 'in', width = 3, height = 0.6)
ggsave('./misc_plots/num_ones_colormap_full.png', units = 'in', width = 3, height = 0.6)

# Full range based on restraint value
ggplot(df, aes(x = num_ones - 60, y = 1, fill = as.factor(num_ones))) + 
  geom_tile() + 
  geom_vline(aes(xintercept = -10), linetype = 'dashed', color = '#000000') +
  scale_fill_manual(values = df$color) +
  scale_x_continuous(expand = c(0,0), breaks = seq(-60,40,10)) +
  scale_y_continuous(expand = c(0,0)) +
  xlab('Restraint value') +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(legend.position = 'none') 
ggsave('./misc_plots/restraint_value_colormap_full.pdf', units = 'in', width = 3, height = 0.6)
ggsave('./misc_plots/restraint_value_colormap_full.png', units = 'in', width = 3, height = 0.6)
 
# Trimmed range based on restraint value
df_plot = df[df$num_ones >= 40 & df$num_ones <= 80,]
ggplot(df_plot, aes(x = num_ones - 60, y = 1, fill = as.factor(num_ones))) + 
  geom_tile() + 
  geom_vline(aes(xintercept = -10), linetype = 'dashed', color = '#000000') +
  scale_fill_manual(values = df_plot$color) +
  scale_x_continuous(expand = c(0,0), breaks = seq(-60,40,10)) +
  scale_y_continuous(expand = c(0,0)) +
  xlab('Restraint value') +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(legend.position = 'none') 
ggsave('./misc_plots/restraint_value_colormap_trimmed.pdf', units = 'in', width = 3, height = 0.6)
ggsave('./misc_plots/restraint_value_colormap_trimmed.png', units = 'in', width = 3, height = 0.6)
 


# Original c++ code

#if(num_ones < threshold){ // Unrestrained
#  color_fill = emp::ColorRGB(
#    255 - (threshold - num_ones) * (255 / threshold),
#    0,
#    0 + (threshold - num_ones) *   (100 / threshold));
#}
#else{ // Restrained
#  color_fill = emp::ColorRGB(
#    255 - (num_ones - threshold) * (255 / (100 - threshold)),
#    255 - (num_ones - threshold) * (255 / (100 - threshold)),
#    255 - (num_ones - threshold) * (255 / (100 - threshold)));
#}