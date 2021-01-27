rm(list = ls())

library(ggplot2)
library(dplyr)


setwd("~/Research/code/SpatialRestraint/scripts/")
cur_size = 8

data = read.csv(paste0('raw_data/', cur_size, '/1_cellGen.dat'))
data_grouped = dplyr::group_by(data, pct_full, cell_generation)
data_summary = dplyr::summarize(data_grouped, num_cells_mean = mean(num_cells))

ggplot(data_summary, aes(x = cell_generation, y = num_cells_mean)) + 
  geom_line() + 
  facet_grid(rows = vars(pct_full)) +
  ggtitle(paste0('Multicell size ', cur_size)) + 
  ggsave(paste0('./plots/', cur_size, '/cell_gens.png'))


ggplot(data, aes(x = cell_generation, y = num_cells, group = X.mc_id)) + 
  geom_line(alpha = 0.25) + 
  facet_grid(rows = vars(pct_full)) +
  ggtitle(paste0('Multicell size ', cur_size)) +
  ggsave(paste0('./plots/', cur_size, '/cell_gens_all_lines.png'))

for(cur_id in 1:5) {
  # Create new data frame to store cell generation median for each pct_full 
  data_medians = data.frame(data = matrix(ncol = 2, nrow = 0))
  # Subset data for the specified mc
  data_mc = data[data$X.mc_id == cur_id,]
  # For each value of pct_full in this mc
  for(pct in unique(data_mc$pct_full)){
    # Create a vector of cell generations (one-to-one, not we have x at y gen)
    gen_vec = c()
    # Subset mc data for this value of pct_full
    data_mc_pct = data_mc[data_mc$pct_full == pct,]
    # For each value of cell_generation present
    for(gen in min(data_mc_pct$cell_generation):max(data_mc_pct$cell_generation)){
      # Subset the data again, for that specific cell_generation
      data_mc_pct_gen = data_mc_pct[data_mc_pct$cell_generation == gen,]
      # Figure out how many cells had that generation value, let's say n cells
      # Add n copies of that generation value to gen_vec
      if(nrow(data_mc_pct_gen) > 0){
        gen_vec = c(gen_vec, rep(gen, data_mc_pct_gen$num_cells[1]))
      }
    }
    # Add new row to the data frame for this value of pct_full
    # Taking the median of the cell generations
    data_medians[nrow(data_medians) + 1, ] = c(pct, median(gen_vec))
  }
  # Set the column names to make our lives easier!
  colnames(data_medians) = c('pct_full', 'cell_generation_median')
  
  # Plot!
  ggplot(data[data$X.mc_id == cur_id,], aes(x = cell_generation, y = num_cells)) + 
    geom_col() +
    geom_vline(data = data_medians, aes(xintercept=cell_generation_median), color = 'red') +
    geom_label(data = data_medians, aes(x = cell_generation_median, y = 1, label=as.character(cell_generation_median))) +
    facet_grid(rows = vars(pct_full)) +
    ggtitle(paste0('Multicell size ', cur_size, ' Multicell #', cur_id)) +
    ggsave(paste0('./plots/', cur_size, '/cell_gens_individual_', cur_id, '.png'))
}


# Start toward cumulative plot
rm(data_grouped)
rm(data_summary)
data_grouped = dplyr::group_by(data, pct_full, cell_generation)
data_sum = dplyr::summarize(data_grouped, num_cells = sum(num_cells))


# Create new data frame to store cell generation median for each pct_full 
data_medians = data.frame(data = matrix(ncol = 2, nrow = 0))
# Subset data for the specified mc
data_mc = data_sum
# For each value of pct_full in this mc
for(pct in unique(data_mc$pct_full)){
  # Create a vector of cell generations (one-to-one, not we have x at y gen)
  gen_vec = c()
  # Subset mc data for this value of pct_full
  data_mc_pct = data_mc[data_mc$pct_full == pct,]
  # For each value of cell_generation present
  for(gen in min(data_mc_pct$cell_generation):max(data_mc_pct$cell_generation)){
    # Subset the data again, for that specific cell_generation
    if(nrow(data_mc_pct[data_mc_pct$cell_generation == gen,]) > 0){
      data_mc_pct_gen = data_mc_pct[data_mc_pct$cell_generation == gen,]
      # Figure out how many cells had that generation value, let's say n cells
      # Add n copies of that generation value to gen_vec
      gen_vec = c(gen_vec, rep(gen, data_mc_pct_gen$num_cells[1]))
    }
  }
  # Add new row to the data frame for this value of pct_full
  # Taking the median of the cell generations
  data_medians[nrow(data_medians) + 1, ] = c(pct, median(gen_vec))
}
# Set the column names to make our lives easier!
colnames(data_medians) = c('pct_full', 'cell_generation_median')

# Plot!
ggplot(data_sum, aes(x = cell_generation, y = num_cells)) + 
  geom_col() +
  geom_vline(data = data_medians, aes(xintercept=cell_generation_median), color = 'red') +
  geom_label(data = data_medians, aes(x = cell_generation_median, y = 5, label=as.character(cell_generation_median))) +
  facet_grid(rows = vars(pct_full)) +
  ggtitle(paste0('Multicell size ', cur_size)) +
  ggsave(paste0('./plots/', cur_size, '/cell_count_sum.png'))

