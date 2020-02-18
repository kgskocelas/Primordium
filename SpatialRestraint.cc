#include <iostream>
#include <fstream>
#include "../../Empirical/source/base/array.h"
#include "../../Empirical/source/base/vector.h"
#include "../../Empirical/source/tools/Random.h"

// Informatin needed to configure a run.
struct Config {
  size_t threshold = 16;
  bool restrain = false;
  size_t neighbors = 8;
};

// Test a single multicell for how long it takes to replicate.
// Multicell is WIDTH * HEIGHT in size
// threshold indicates resources needed to replicate.
// restrain indicates if cells should hold back.

template <size_t WIDTH, size_t HEIGHT>
struct World {
  static constexpr size_t SIZE = WIDTH * HEIGHT;
    
  static constexpr size_t ToPos(size_t x, size_t y) { return x + y * WIDTH; }
  static constexpr size_t ToX(size_t pos) { return pos % WIDTH; }
  static constexpr size_t ToY(size_t pos) { return pos / WIDTH; }

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
  // Thus 0-1 is a 1D size 2 neighborhood; 0-3 are a 2D size-4 neighborhood; 0-7 is a 2D size 8 neighborhood.
  // (0-5 behaves like a hex-map) Larger assumes popoulation size and returns the full set.

  static size_t RandomNeighbor(emp::Random & random, size_t pos, size_t neighbors=8)
  {
    if (neighbors > 8) {
      emp_assert(neighbors == SIZE, neighbors, SIZE);
      return random.GetUInt(neighbors);
    }

    const size_t x = ToX(pos);
    const size_t y = ToY(pos);
    size_t next_x = (size_t) -1;
    size_t next_y = (size_t) -1;

    while (next_x >= WIDTH || next_y >= HEIGHT) {
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

    emp_assert(ToPos(next_x, next_y) < SIZE, next_x, next_y, WIDTH, HEIGHT);

    return ToPos(next_x, next_y);
  }

  static void Print(const emp::array<size_t, SIZE> & ms) {
    size_t pos = 0;
    for (size_t y = 0; y < HEIGHT; y++) {
      for (size_t x = 0; x < WIDTH; x++) {
	std::cout << " " << ToChar(ms[pos++]);
      }
      std::cout << std::endl;
    }
  }
  
  static size_t TestMulticell(emp::Random & random, size_t threshold, bool restrain, size_t neighbors=8) {
    // Setup initial multicell to be empty; keep count of resources in each cell.
    emp::array<size_t, SIZE> orgs;
    for (size_t i = 0; i < SIZE; i++) orgs[i] = 0;
    size_t time = 0;
    
    // Inject a cell in the middle.
    constexpr size_t start_pos = ToPos(WIDTH/2, HEIGHT/2);
    orgs[start_pos] = 1; // Live cells start at one resource.
    size_t num_orgs = 1;    

    // Loop through updates until cell is full.
    while (num_orgs < SIZE) {
      // std::cout << "Update: " << time << "; num_orgs: " << num_orgs << std::endl;
      // if (time % 100 == 0) Print(orgs);

      time++;

      // Loop through cells, probabilistically incrementing resources.
      for (size_t pos = 0; pos < SIZE; pos++) {
        // See if this cell gets a unit of resource.
        if (orgs[pos] && random.P(0.5)) {
          // If so, give resource and see if it replicates.
          if (++orgs[pos] == threshold) {
            orgs[pos] = 1;

            size_t next_pos = RandomNeighbor(random, pos, neighbors);
            size_t & next_org = orgs[next_pos];

            // If the target is empty, put a new organism there.
            if (!next_org) {
              next_org = 1;
              num_orgs++;
            }

            // Otherwise if we don't restrain, reset organism there.
            else if (!restrain) {
              next_org = 1;
            }
          }
        }
      }
    }

    return time;
  }
};

// General case
template <size_t...> struct WorldSet;

// Base case
template <>
struct WorldSet<> {
  template <typename... Ts>
  static void Run(Ts &&...) { ; }
};

// Recursive case
template <size_t CUR_SIZE, size_t... WORLD_SIZES>
struct WorldSet<CUR_SIZE, WORLD_SIZES...> {
  static void Run(emp::Random & random,
		  std::ostream & os=std::cout,
		  size_t threshold=20,
      size_t neighbors=8,
		  size_t num_runs=100,
		  bool verbose=false) {
    // Build a world of the correct size.
    World<CUR_SIZE,CUR_SIZE> world;
    
    // Output initial size information, first for no restraint
    os << CUR_SIZE << ", " << CUR_SIZE;
    if (verbose) os << ", 0";
    
    double total_time = 0.0;
    for (size_t i = 0; i < num_runs; i++) {
      size_t cur_time = world.TestMulticell(random, threshold, false, neighbors);
      if (verbose) os << ", " << cur_time;
      total_time += (double) cur_time;
    }
    os << ", " << (total_time / (double) num_runs);
    if (verbose) os << std::endl;
    
    // Then for restraint...
    if (verbose) os << CUR_SIZE << ", " << CUR_SIZE << ", 1";
    
    total_time = 0.0;
    for (size_t i = 0; i < num_runs; i++) {
      size_t cur_time = world.TestMulticell(random, threshold, true, neighbors);
      if (verbose) os << ", " << cur_time;
      total_time += (double) cur_time;
    }
    os << ", " << (total_time / (double) num_runs) << std::endl;
    
    
    // Do the recursive call.
    WorldSet<WORLD_SIZES...>::Run(random, os, threshold, neighbors, num_runs, verbose);
  }
};

int main()
{
  emp::Random random;
  WorldSet<2,4,8,16,32,64,128>::Run(random, std::cout, 20, 8, 100, true);
}
