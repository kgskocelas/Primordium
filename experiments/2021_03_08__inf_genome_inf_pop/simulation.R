rm(list = ls())

# Load libraries
library(ggplot2)
library(dplyr)

# Load data
df_summary = read.csv('./timing/data/timing_data_summary.csv')

# Collect fitness data in an easier-to-access format
size_vec = c(16, 32, 64, 128, 256, 512, 1024)
df_fitness = data.frame(data = matrix(nrow = 0, ncol = 3))
colnames(df_fitness) = c('mc_size', 'num_ones', 'fitness')
for(mc_size in size_vec){
  size_mask = df_summary$cells_side == mc_size
  for(ones in unique(df_summary$ancestor_1s)){
    df_fitness[nrow(df_fitness) + 1,] = c(mc_size, ones, df_summary[size_mask & df_summary$ancestor_1s == ones,]$mean_rep_rate)
  }
}

# Create directory to save data
if(!dir.exists('./simulation')){
  dir.create('./simulation')
}
if(!dir.exists('./simulation/data')){
  dir.create('./simulation/data')
}
# Create directory for plots
if(!dir.exists('./simulation/plots')){
  dir.create('./simulation/plots')
}
if(!dir.exists('./simulation/plots/sim')){
  dir.create('./simulation/plots/sim')
}

# Global variables for simulation
num_generations = 10000
gen_write_step = 500
mut_rate = 0.02
mut_pos = mut_rate * 0.4
mut_neg = mut_rate * 0.6

# Iterate through each organism size
for(size_idx in 1:length(size_vec)){
  size = size_vec[size_idx]
  genome_length = 650
  # Create data frame to collect simulated data
  df_density = data.frame(data = matrix(nrow = (genome_length + 1), ncol = 4))
  colnames(df_density) = c('mc_size', 'num_ones', 'generation', 'frac_of_pop')
  gen_offset = genome_length + 1
  fit_size_mask = df_fitness$mc_size == size
  fit_vec = df_fitness[fit_size_mask,]$fitness 
  cat('Size: ', size, '\n')
  # Create vectors to hold previous and current population fractions
  prev_fracs = rep(0, genome_length + 1)
  prev_fracs[100 + 1] = 1
  cur_fracs = rep(0, genome_length + 1)
  for(num_ones in 0:genome_length){
    df_density[num_ones + 1,] = c(size, num_ones - 100, 0, prev_fracs[num_ones+ 1])
  }
  # Iterate through each 
  for(gen in 1:num_generations){
    for(num_ones in 0:genome_length){
      if(num_ones == 0){ # Lower edge case
        mut_rate_from_next = mut_neg
        mut_rate_to_next = mut_pos
        cur_fracs[num_ones + 1] = 
            (1 - mut_rate_to_next) * (fit_vec[num_ones + 1]) * prev_fracs[num_ones + 1] + 
            (fit_vec[num_ones + 2]) * prev_fracs[num_ones + 2] * mut_rate_from_next 
      }else if(num_ones == genome_length){ # Upper edge case
        mut_rate_from_prev = mut_pos
        mut_rate_to_prev = mut_neg
        cur_fracs[num_ones + 1] = 
            (1  - mut_rate_to_prev) * (fit_vec[num_ones + 1]) * prev_fracs[num_ones + 1] + 
            (fit_vec[num_ones]) * prev_fracs[num_ones] * mut_rate_from_next 
      }else{ # Base case
        mut_rate_from_next = mut_neg
        mut_rate_to_next =   mut_pos
        mut_rate_from_prev = mut_pos
        mut_rate_to_prev =   mut_neg
        cur_fracs[num_ones + 1] = 
            (1 - mut_rate_to_next - mut_rate_to_prev) * (fit_vec[num_ones + 1]) *prev_fracs[num_ones + 1] + 
            (fit_vec[num_ones + 2]) * prev_fracs[num_ones + 2] * mut_rate_from_next + 
            (fit_vec[num_ones]) * prev_fracs[num_ones] * mut_rate_from_prev
      }
    }
    # Normalize fraction to be a perfect 1.0
    cur_fracs = cur_fracs / sum(cur_fracs)
    # Record ever k generations 
    if(gen %% gen_write_step == 0){
      cat(gen, ' ')
      for(num_ones in 0:genome_length){
        df_density[nrow(df_density) + 1,] = c(size, num_ones - 100, gen, cur_fracs[num_ones+ 1])
      }
    }
    # Save current population fractions for next step
    prev_fracs = cur_fracs
  }
  cat('\n')
 
  # Plot simulated data
  mask = !is.na(df_density$mc_size) & df_density$mc_size == size
  ggplot(df_density[mask,], aes(x = num_ones, y = generation, fill = frac_of_pop)) + 
    geom_raster() + 
    scale_y_reverse(expand = c(0,0)) + 
    scale_x_continuous(expand = c(0,0)) 
  ggsave(paste0('./simulation/plots/sim/density_over_time__size_', size, '.png'), units = 'in', width = 6, height = 6)
  
  # Save simulated data
  write.csv(df_density, paste0('simulation/data/size_', size, '.csv'))
  
  mask = !is.na(df_density$mc_size) & df_density$mc_size == size & df_density$generation == num_generations
  cat('Mean at last generation: ', weighted.mean(-100:550, df_density[mask,]$frac_of_pop) ,'\n')
}
