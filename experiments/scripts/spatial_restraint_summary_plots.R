rm(list = ls())
library(ggplot2)

#### Begin Configuration
args = commandArgs(trailingOnly=T)
if(length(args) != 2){
    print('Error! Must pass exactly 6 command line arguments:')
    print('1. Input filename')
    print('2. Output directory')
    q()
}
# File to plot
input_filename = args[1]
# Where to save plots 
output_dir = args[2]
# Appended to the end of filenames for both loading and saving to make organization easier
#suffix = 'example'
# Dimensions of output in inches
width  = 8
height = 6
# Do we output pdf/png?
output_pdf = TRUE
output_png = FALSE
# Do we want to create boxplots for each ancestor starting point?
multiple_starts = FALSE
# What costs (in terms of updates per unrestrained cell at replication) should we plot?
cost_vec = c(0, 100)
# Should we draw each ancestral number of ones (50, 75, 100, etc) with different line style?
use_linetype_ones = FALSE
# Save plots with 95% confidence interval ribbons? (in addition)
use_ribbon = TRUE
# If true, also save plots with only the selected (see next var) multicell sizes
use_selected_sizes = TRUE
# If use_selected_sizes is true, these are the multicell sizes to plot (thinned to cur clutter)
selected_sizes_vec = c(512, 128, 32)
#### End Configuration


#### Load in the summary data
print('Loading data, this may take a moment.')
data = read.csv(input_filename)#paste0('./intermediate_data/data_', suffix, '_summary.csv'))
print('Done loading data!')

#### Begin plots
for(cur_cost in cost_vec){
  data_plot = data[data$COST == cur_cost,]
  if(nrow(data_plot) < 1){
    print(paste0('No rows for cost = ', cur_cost, '. Skipping!'))
    next
  }
  # Create directory for this cost level (historical contingency)
  dir.create(output_dir)
  dir_cost = output_dir
  
  
  # Plot the grand mean for each multicell size (and optionally ancestor point)
  ggp = ggplot(data_plot, 
    aes(x = generation, y = grand_mean_ones, color = as.factor(MCSIZE)))
  if(use_linetype_ones){
    ggp = ggp + geom_line(aes(linetype = as.factor(ONES)))
  }else{
    ggp = ggp + geom_line()
  }
  ggp = ggp +
    geom_line(aes(linetype = as.factor(ONES))) + 
    xlab('Generation') + 
    ylab('Number of ones') +
    labs(color = 'Cells per side') + 
    ggtitle('No extra cost for unrestrained cells')
  if(output_pdf){
    ggp + ggsave(paste0(dir_cost, '/grand_mean.pdf'), 
      units = 'in', width = width, height = height)
  }
  if(output_png){
    ggp + ggsave(paste0(dir_cost, '/grand_mean.png'), 
      units = 'in', width = width, height = height)
  }


  # Plot the grand mean for each multicell size (and optionally ancestor point) with ribbons
  if(use_ribbon){
    ggp = ggplot(data_plot, 
      aes(x = generation, y = grand_mean_ones, color = as.factor(MCSIZE)))
    ggp = ggp + geom_ribbon(
      aes(ymin = grand_mean_ones - offset_95, ymax = grand_mean_ones + offset_95), alpha = 0.2)
    if(use_linetype_ones){
      ggp = ggp + geom_line(aes(linetype = as.factor(ONES)))
    }else{
      ggp = ggp + geom_line()
    }
    ggp = ggp +
      geom_line(aes(linetype = as.factor(ONES))) + 
      xlab('Generation') + 
      ylab('Number of ones') +
      labs(color = 'Cells per side') + 
      ggtitle('No extra cost for unrestrained cells')
    if(output_pdf){
      ggp + ggsave(paste0(dir_cost, '/grand_mean_ribbon.pdf'), 
        units = 'in', width = width, height = height)
    }
    if(output_png){
      ggp + ggsave(paste0(dir_cost, '/grand_mean_ribbon.png'), 
        units = 'in', width = width, height = height)
    }
  }
 
  # Should we also create some thinner plots? 
  if(use_selected_sizes){
    data_plot_selected = data_plot[data_plot$MCSIZE %in% selected_sizes_vec,]
    # Plot the grand mean for each selected multicell size (and optionally ancestor point)
    ggp = ggplot(data_plot_selected, 
      aes(x = generation, y = grand_mean_ones, color = as.factor(MCSIZE)))
    if(use_linetype_ones){
      ggp = ggp + geom_line(aes(linetype = as.factor(ONES)))
    }else{
      ggp = ggp + geom_line()
    }
    ggp = ggp +
      geom_line(aes(linetype = as.factor(ONES))) + 
      xlab('Generation') + 
      ylab('Number of ones') +
      labs(color = 'Cells per side') + 
      ggtitle('No extra cost for unrestrained cells')
    if(output_pdf){
      ggp + ggsave(paste0(dir_cost, '/grand_mean_selected.pdf'), 
        units = 'in', width = width, height = height)
    }
    if(output_png){
      ggp + ggsave(paste0(dir_cost, '/grand_mean_selected.png'), 
        units = 'in', width = width, height = height)
    }


    # Plot the grand mean for each selected multicell size (and optionally ancestor point) w/ ribbons
    if(use_ribbon){
      ggp = ggplot(data_plot_selected, 
        aes(x = generation, y = grand_mean_ones, color = as.factor(MCSIZE)))
      ggp = ggp + geom_ribbon(
        aes(ymin = grand_mean_ones - offset_95, ymax = grand_mean_ones + offset_95), alpha = 0.2)
      if(use_linetype_ones){
        ggp = ggp + geom_line(aes(linetype = as.factor(ONES)))
      }else{
        ggp = ggp + geom_line()
      }
      ggp = ggp +
        geom_line(aes(linetype = as.factor(ONES))) + 
        xlab('Generation') + 
        ylab('Number of ones') +
        labs(color = 'Cells per side') + 
        ggtitle('No extra cost for unrestrained cells')
      if(output_pdf){
        ggp + ggsave(paste0(dir_cost, '/grand_mean_selected_ribbon.pdf'), 
          units = 'in', width = width, height = height)
      }
      if(output_png){
        ggp + ggsave(paste0(dir_cost, '/grand_mean_selected_ribbon.png'), 
          units = 'in', width = width, height = height)
      }
    }
  }
}
