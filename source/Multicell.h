/**
 *  @copyright Copyright (C) Michigan State University, MIT Software license; see doc/LICENSE.md
 *  @date 2020.
 *
 *  @file  Multicell.h
 *  @brief Manages analysis of a single multicell
 *  @note Status: BETA
 */

#ifndef MULTICELL_H
#define MULTICELL_H

#include "base/vector.h"
#include "tools/Random.h"
#include "tools/stats.h"
#include "tools/TimeQueue.h"

/// Information about a single cell.
struct Cell {
  size_t id;
  double repro_time = 0.0;  ///< When will this cell replicate?
  size_t num_ones = 0;      ///< How many ones in genome?

  bool operator==(const Cell & _in) const { return id == _in.id; }
  bool operator!=(const Cell & _in) const { return id != _in.id; }
  bool operator<(const Cell & _in) const {
    if (repro_time == _in.repro_time) return (id < _in.id);
    return (repro_time < _in.repro_time);
  }
};

/// Results from a single run.
struct RunResults {
  double run_time;                  ///< What was the replication time of this group?
  emp::vector<double> cell_counts;  ///< How many cells have each bit count?
  double extra_cost;                ///< Extra cost due to unrestrained cells.

  RunResults() : run_time(0.0), cell_counts(0) { ; }
  RunResults(const size_t num_bits) : run_time(0.0), cell_counts(num_bits+1, 0.0) { ; }
  RunResults(const RunResults &) = default;
  RunResults(RunResults &&) = default;

  RunResults & operator=(const RunResults &) = default;
  RunResults & operator=(RunResults &&) = default;

  RunResults & operator+=(const RunResults & in) {
    emp_assert(cell_counts.size() == in.cell_counts.size());
    run_time += in.run_time;
    for (size_t i=0; i < cell_counts.size(); i++) cell_counts[i] += in.cell_counts[i];
    extra_cost += extra_cost;
    return *this;
  }

  RunResults & operator/=(const double denom) {
    emp_assert(denom != 0.0);
    run_time /= denom;
    for (double & x : cell_counts) x /= denom;
    extra_cost /= denom;
    return *this;
  }

  /// Count the total number of cells represented.
  double CountCells() const { return emp::Sum(cell_counts); }

  /// Count the number of cells that exhibit restrained behavior.
  double CountRestrained(size_t threshold) const {
    double total = 0.0;
    for (size_t i = threshold; i < cell_counts.size(); i++) total += cell_counts[i];
    return total;
  }

  /// Count the number of cells that DO NOT exhibit restrained behavior.
  double CountUnrestrained(size_t threshold) const {
    double total = 0.0;
    for (size_t i = 0; i < threshold; i++) total += cell_counts[i];
    return total;
  }

  /// Caclulate the full replication time.
  double GetReproTime() const { return run_time + extra_cost; }
};

/// A single "multicell" organism.

struct Multicell {
  emp::Random & random;

  emp::vector<Cell> cells;   ///< All cells in this multicell
  emp::vector<char> is_full; ///< Is the local neighborhood full?
  size_t num_cells = 0;      ///< How many cells are currnetly in the multicell?
  size_t mask_side = 31;     ///< Bit mask for a side (for id -> x pos)
  size_t log2_side = 5;      ///< Log base 2 of the number of cells on a side (for id -> y pos).

  emp::TimeQueue<size_t> cell_queue; ///< Cells waiting to replicate.

  size_t time_range = 50.0;  ///< Replication takes 100.0 + a random value up to time_range.
  size_t neighbors = 8;      ///< Num neighbors in grid for offspring (0=well mixed; 4,6,8 => 2D)
  size_t cells_side = 32;    ///< How many cells are on a side of the (square) multi-cell?
  size_t genome_size = 10;   ///< How many bits in genome?
  size_t restrain = 5;       ///< How many ones in bit sequence for restraint?
  size_t start_1s = 5;       ///< How many ones in the starting cell?
  double mut_prob = 0.0;     ///< Probability of an offspring being mutated.
  double unrestrained_cost = 0.0; ///< Extra cost for each unrestrained cell when full.
  bool one_check = false;    ///< Should restrained check only one cell to find empty?

  Multicell(emp::Random & _random) : random(_random), cell_queue(100.0) { }

  size_t GetSize() const { return cells_side * cells_side; }

  size_t ToPos(size_t x, size_t y) const { return x + y * cells_side; }
  size_t ToX(size_t pos) const { return pos & mask_side; }
  size_t ToY(size_t pos) const { return pos >> log2_side; }

  size_t MiddlePos() const { return ToPos(cells_side/2, cells_side/2); }

  // Convert a resource count to a character.
  static constexpr char ToChar(size_t count) {
    if (count < 10) return '0' + count;
    if (count < 36) return 'a' + (count-10);
    if (count < 62) return 'A' + (count-36);
    return '+';
  }

  // Neighborhood layout:
  //  7 2 4
  //  0 * 1
  //  5 3 6
  //
  // Thus 0-1 is a 1D size 2 neighborhood; 0-3 are a 2D size-4 neighborhood; 0-7 is a 2D size 8
  // neighborhood.  (0-5 behaves like a hex-map) Larger assumes popoulation size and returns the
  // full set.

  size_t RandomNeighbor(size_t pos)
  {
    if (neighbors == 0 || neighbors > 8) {
      return random.GetUInt(GetSize());
    }

    const size_t x = ToX(pos);
    const size_t y = ToY(pos);
    size_t next_x = (size_t) -1;
    size_t next_y = (size_t) -1;

    while (next_x >= cells_side || next_y >= cells_side) {
      const size_t dir = random.GetUInt(neighbors);  // Direction for offspring.
      switch (dir) {
      case 0: case 5: case 7: next_x = x-1; break;
      case 1: case 4: case 6: next_x = x+1; break;
      default: next_x = x;
      };
      switch (dir) {
      case 2: case 4: case 7: next_y = y-1; break;
      case 3: case 5: case 6: next_y = y+1; break;
      default: next_y = y;
      };
    }

    emp_assert(ToPos(next_x, next_y) < GetSize(), next_x, next_y, cells_side);

    return ToPos(next_x, next_y);
  }

  size_t EmptyNeighbor(size_t pos)
  {
    if (is_full[pos] == 1) return (size_t) -1;

    // If well mixed, keep searching until we find a value.
    if (neighbors == 0 || neighbors > 8) {
      size_t id = random.GetUInt(GetSize());
      while (cells[id].repro_time != 0.0) id = random.GetUInt(GetSize());
      return id;
    }

    // Use a static vector so that we don't keep reallocating.
    static emp::vector<size_t> neighbor_ids;
    neighbor_ids.resize(0);

    const size_t x = ToX(pos);
    const size_t y = ToY(pos);
    size_t next_x = (size_t) -1;
    size_t next_y = (size_t) -1;

    for (size_t dir = 0; dir < neighbors; dir++) {
      switch (dir) {
      case 0: case 5: case 7: next_x = x-1; break;
      case 1: case 4: case 6: next_x = x+1; break;
      default: next_x = x;
      };
      switch (dir) {
      case 2: case 4: case 7: next_y = y-1; break;
      case 3: case 5: case 6: next_y = y+1; break;
      default: next_y = y;
      };

      // If this position is on the grid AND empty
      if (next_x < cells_side && next_y < cells_side) {
        const size_t next_pos = ToPos(next_x, next_y);
        if (cells[next_pos].repro_time == 0.0) neighbor_ids.push_back(next_pos);
      }
    }

    if (neighbor_ids.size() == 0) {
      is_full[pos] = 1;
      return (size_t) -1;
    }

    return neighbor_ids[random.GetUInt(neighbor_ids.size())];
  }

  // Print current one-counts in population.
  void Print() {
    emp_assert(cells.size() == GetSize());
    size_t pos = 0;
    for (size_t y = 0; y < cells_side; y++) {
      for (size_t x = 0; x < cells_side; x++) {
        if (cells[pos].repro_time == 0.0) std::cout << " -";
      	else std::cout << " " << ToChar(cells[pos].num_ones);
        pos++;
      }
      std::cout << std::endl;
    }
  }

  void SetupCell(Cell & cell) {
    cell.repro_time = cell_queue.GetTime() + 100.0 + random.GetDouble(time_range);
    cell_queue.Insert(cell.id, cell.repro_time);
  }

  void InjectCell(size_t pos, size_t num_ones) {
    Cell & inject_cell = cells[pos];                 // Find cell at inject position.
    if (inject_cell.repro_time == 0.0) num_cells++;  // If cell was empty, mark increase.
    inject_cell.num_ones = num_ones;                 // Initialize injection ones.
    SetupCell(inject_cell);                          // Do any extra setup for this cell.
  }

  void InjectCell(size_t pos) { InjectCell(pos, start_1s); }

  // Setup the new offspring, possibly with mutations.
  void DoBirth(Cell & offspring, const Cell & parent, bool do_mutations=true) {
    if (offspring.repro_time == 0.0) num_cells++;  // If offspring was empty, this is a new cell.
    offspring.num_ones = parent.num_ones;
    if (do_mutations && random.P(mut_prob)) {
      double prob1 = ((double) offspring.num_ones) / (double) genome_size;
      if (random.P(prob1)) offspring.num_ones--;
      else offspring.num_ones++;
    }

    SetupCell(offspring);       // Launch cell in the population.
    is_full[offspring.id] = 0;  // Mark local region as NOT FULL.
  }

  /// Once we have current settings locked in, reset all non-setting values appropriately.
  void SetupConfig() {
    // Setup initial multicell to be empty; keep count of resources in each cell.
    cells.resize(0);           // Clear out any current cells.
    cells.resize(GetSize());   // Put in new cells initialized to 0.
    for (size_t id = 0; id < cells.size(); id++) {
      cells[id].id = id;
      cells[id].repro_time = 0.0;
    }
    is_full.resize(0);
    is_full.resize(GetSize(), 0);
    cell_queue.Reset();
    num_cells = 0;

    if (emp::count_bits(cells_side) != 1) {
      std::cerr << "\nERROR: Cannot have " << cells_side << "cells on a side; must be a power of 2!\n";
      exit(1);
    }
    mask_side = cells_side - 1;
    log2_side = emp::count_bits(mask_side);
  }

  /// Run the multicell until it is full.
  RunResults Run(bool print_trace=false, std::ostream & os=std::cout) {
    size_t last_count = 0;                   // Track cells from last time (for traces)

    while (num_cells < cells.size()) {
      emp_assert(cell_queue.GetSize() > 0);

      Cell & parent = cells[cell_queue.Next()];

      // If this cell has been updated since being bufferred, skip it.
      if (parent.repro_time != cell_queue.GetTime()) continue;

      // Neighborhood is only marked full for restrained orgs; if so, fail divide.
      if (is_full[parent.id]) continue; 

      size_t next_id = RandomNeighbor(parent.id); // Find the placement of the offspring.
      Cell & next_cell = cells[next_id];

      // If the target is empty or we don't restrain, put a new cell there.
      if (next_cell.repro_time == 0.0 || parent.num_ones < restrain) {
        DoBirth(next_cell, parent);
      }

      // Otherwise it is restrained and not empty; unless limited  to one, keep looking!
      else if (!one_check) {
        next_id = EmptyNeighbor(parent.id);
        if (next_id != (size_t) -1) DoBirth(cells[next_id], parent);
      }

      SetupCell(parent);  // Reset parent for its next replication.

      // If we are tracing, output data.
      if (print_trace && last_count != num_cells) {
        last_count = num_cells;
        os << "\nTime: " << cell_queue.GetTime()
           << "  Cells: " << last_count
           << "\n";
        Print();
      }
    }

    // Setup the results and return them.
    RunResults results(genome_size);
    results.run_time = cell_queue.GetTime();
    size_t unrestrained_count = 0;
    for (const auto & cell : cells) {
      if (cell.num_ones < restrain) unrestrained_count++;
      results.cell_counts[cell.num_ones] += 1.0;
    }
    results.extra_cost = unrestrained_count * unrestrained_cost;

    return results;
  }
};

#endif
