rm(list = ls())
library(ggplot2)

#### Begin Configuration
args = commandArgs(trailingOnly=T)
if(length(args) != 3){
    print('Error! Must pass exactly 6 command line arguments:')
    print('1. Input filename')
    print('2. Output directory')
    print('3. Max generation')
    q()
}
# File to plot
input_filename = args[1]
# Where to save plots 
output_dir = args[2]
# Max generation
max_gen = as.numeric(args[3])
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
cost_vec = c(0)
#### End Configuration

#### Load in the cleaned data
print('Loading data, this may take a moment.')
data = read.csv(input_filename)#paste0('./intermediate_data/data_', suffix, '.csv'))
print('Done loading data!')
print(paste0('Saving plots to ', output_dir))

#### Begin plots
for(cur_cost in cost_vec){
  # Prepare data for boxplots
  data_plot = data[data$generation == max_gen & data$COST == cur_cost,]
  if(nrow(data_plot) < 1){
    print(paste0('No rows for cost = ', cur_cost, '. Skipping!'))
    next
  }
  data_plot$MCSIZE = as.factor(as.numeric(data_plot$MCSIZE))
  
  # Create directory for this cost level (historical contingency)
  dir.create(output_dir)
  dir_cost = output_dir 
  
    
  # Boxplot of average ones by multicell size
  ggp = ggplot(data_plot, aes(x = as.factor(MCSIZE), y = ave_ones)) + 
    geom_boxplot() + 
    xlab('Cells per multicells side') + 
    ylab('Average number of ones in genome')
  if(output_pdf){
    ggp + ggsave(paste0(dir_cost, '/boxplot_ones.pdf'), 
           units = 'in', width = width, height = height)
  }
  if(output_png){
    ggp + ggsave(paste0(dir_cost, '/boxplot_ones.png'), 
           units = 'in', width = width, height = height)
  }
  rm(ggp)


  # Boxplot of average reproduction time by multicell size
  ggp = ggplot(data_plot, aes(x = as.factor(MCSIZE), y = ave_repro_time, fill =as.factor(MCSIZE))) +
    geom_boxplot() +
    xlab('Cells per multicells side') + 
    ylab('Multicell reproduction time')
  if(output_pdf){
    ggp + ggsave(paste0(dir_cost, '/boxplot_repro_time.pdf'), 
           units = 'in', width = width, height = height)
  }
  if(output_png){
    ggp + ggsave(paste0(dir_cost, '/boxplot_repro_time.png'), 
           units = 'in', width = width, height = height)
  }
  rm(ggp)


  if(multiple_starts){
    # Boxplot of average ones by multicell size and ancestor starting ones
    ggp = ggplot(data_plot, 
            aes(
                  x = as.factor(MCSIZE), 
                  y = ave_ones, 
                  fill = as.factor(ONES)), 
            group = as.factor(ONES)) + 
      geom_boxplot(position = position_dodge()) + 
      xlab('Cells per multicells side') + 
      ylab('Average number of ones in genome') +
      labs(fill = 'Starting ones')
    if(output_pdf){
      ggp + ggsave(paste0(dir_cost, '/boxplot_ones_starts.pdf'), 
             units = 'in', width = width, height = height)
    }
    if(output_png){
      ggp + ggsave(paste0(dir_cost, '/boxplot_ones_starts.png'), 
             units = 'in', width = width, height = height)
    }
    rm(ggp)
    

    # Boxplot of multicell reproduction time by multicell size and ancestor starting ones
    ggp = ggplot(data_plot, 
            aes(
                  x = as.factor(MCSIZE), 
                  y = ave_repro_time, 
                  fill = as.factor(ONES)), 
            group = as.factor(ONES)) + 
      geom_boxplot(position = position_dodge()) + 
      xlab('Cells per multicells side') + 
      ylab('Multicell reproduction time') +
      labs(fill = 'Starting ones')
    if(output_pdf){
      ggp + ggsave(paste0(dir_cost, '/boxplot_repro_time_starts.pdf'), 
             units = 'in', width = width, height = height)
    }
    if(output_png){
      ggp + ggsave(paste0(dir_cost, '/boxplot_repro_time_starts.png'), 
             units = 'in', width = width, height = height)
    }
    rm(ggp)
  }

  # Violin plot of average ones by multicell size and ancestor starting ones
  ggp = ggplot(data_plot, 
          aes(
                x = as.factor(MCSIZE), 
                y = ave_ones, 
                fill = as.factor(ONES)), 
          group = as.factor(ONES)) + 
    geom_violin() + 
    xlab('Cells per multicells side') + 
    ylab('Average number of ones in genome') +
    labs(fill = 'Starting ones')
  if(output_pdf){
    ggp + ggsave(paste0(dir_cost, '/violin_ones_starts.pdf'), 
           units = 'in', width = width, height = height)
  }
  if(output_png){
    ggp + ggsave(paste0(dir_cost, '/violin_ones_starts.png'), 
           units = 'in', width = width, height = height)
  }
  rm(ggp)
  

  # Violin plot of multicell reproduction time by multicell size and ancestor starting ones
  ggp = ggplot(data_plot, 
          aes(
                x = as.factor(MCSIZE), 
                y = ave_repro_time, 
                fill = as.factor(ONES)), 
          group = as.factor(ONES)) + 
    geom_violin() + 
    xlab('Cells per multicells side') + 
    ylab('Multicell reproduction time') +
    labs(fill = 'Starting ones')
  if(output_pdf){
    ggp + ggsave(paste0(dir_cost, '/violin_repro_time_starts.pdf'), 
           units = 'in', width = width, height = height)
  }
  if(output_png){
    ggp + ggsave(paste0(dir_cost, '/violin_repro_time_starts.png'), 
           units = 'in', width = width, height = height)
  }
  rm(ggp)
}



# There be monsters!

#data$line_id = paste0(data$MCSIZE, '_', data$rep_id, '_', data$ONES)
#ggplot(data[data$MCSIZE == 512 & data$ONES == 100,], aes(x = generation, y = ave_ones, color = as.factor(MCSIZE), group = as.factor(line_id))) + 
#  geom_line(alpha = 0.25) 
#
#ggplot(data[data$MCSIZE == 512,], aes(x = generation, y = ave_ones, color = as.factor(ONES), group = as.factor(line_id))) + 
#  geom_line(alpha = 0.25) +
#  ggsave(paste0('spatial_restraint__512_x_ones__', suffix, '.pdf'), units = 'in', width = 8, height = 6) +
#  ggsave(paste0('spatial_restraint__512_x_ones__', suffix, '.png'), units = 'in', width = 8, height = 6)
#

