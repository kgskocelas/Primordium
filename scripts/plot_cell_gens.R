rm(list = ls())

library(ggplot2)
library(dplyr)

data = read.csv('cell_gens.dat')
data_grouped = dplyr::group_by(data, pct_full, cell_generation)
data_summary = dplyr::summarize(data_grouped, num_cells_mean = mean(num_cells))

ggplot(data_summary, aes(x = cell_generation, y = num_cells_mean)) + 
  geom_line() + 
  facet_grid(rows = vars(pct_full)) +
  ggsave('./scripts/plots/cell_gens.png')


ggplot(data, aes(x = cell_generation, y = num_cells, group = X.mc_id)) + 
  geom_line(alpha = 0.25) + 
  facet_grid(rows = vars(pct_full)) +
  ggsave('./scripts/plots/cell_gens_all_lines.png')
