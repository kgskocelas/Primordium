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

#endif