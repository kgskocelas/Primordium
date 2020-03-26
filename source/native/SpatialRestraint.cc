/**
 *  @copyright Copyright (C) Michigan State University, MIT Software license; see doc/LICENSE.md
 *  @date 2020.
 *
 *  @file  SpatialRestraint.cc
 *  @brief Examining the evolution of cellular restraint in a multi-cell
 *  @note Status: BETA
 *
 *  Cells are either RESTRAINED or UNRESTRAINED
 *  Organisms are "multi-cells" that start as a single cell and automatically replicate when full.
 *  Cells replicate within the multicell when they have collected enough resource.
 *  - Restrained cells replicate into an empty neighbor position (or fail if no empty)
 *  - Unrestrained cells replicate into a random position, regardless of if anything is there.
 */

#include "../SpatialRestraint.h"

int main(int argc, char* argv[])
{
  World world(argc, argv);
  world.Run(std::cout);
}
