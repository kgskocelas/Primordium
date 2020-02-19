#include <iostream>
#include <fstream>
#include "../../Empirical/source/base/array.h"
#include "../../Empirical/source/base/vector.h"
#include "../../Empirical/source/config/command_line.h"
#include "../../Empirical/source/tools/Random.h"
#include "../../Empirical/source/tools/string_utils.h"
#include "../../Empirical/source/tools/vector_utils.h"

// Informatin needed to configure a run.
struct Config {
  bool restrain = false;
  size_t threshold = 16;
  size_t neighbors = 8;

  size_t num_runs = 100;

  std::string GetHeaders() const { return "threshold, restrain, neighbors"; }
  std::string AsCSV() const { return emp::to_string(threshold, ", ", restrain, ", ", neighbors); }
};

// Set of parameters to use.
struct ConfigSet {
  emp::vector<bool> restrain_set;           ///< Should organisms be restrained?
  emp::vector<size_t> threshold_set;        ///< How high a count for replication?
  emp::vector<size_t> neighbor_set;         ///< How many neighbors should we look at; (0=well mixed; 4,6,8 are 2D)

  emp::array<size_t, 3> cur_ids{0,0,0};     ///< Which settings are we going to use next?

  size_t num_runs = 100;                    ///< How many times should we run each configuration?
  size_t max_size = (size_t) -1;            ///< Largest SIDE for a population.

  size_t GetSize() const { return restrain_set.size() * threshold_set.size() * neighbor_set.size(); }

  Config GetConfig() const {
    return Config{ restrain_set[cur_ids[0]], threshold_set[cur_ids[1]], neighbor_set[cur_ids[2]], num_runs };
  }

  bool Next() {
    // First try to toggle restraint.
    cur_ids[0]++;
    if (cur_ids[0] < restrain_set.size()) return true;

    // Next, cycle restrain back to the beginning and try to move to the next threshold.
    cur_ids[0] = 0;
    cur_ids[1]++;
    if (cur_ids[1] < threshold_set.size()) return true;

    // If we're done with thresholds, cycle back and try neighbors.
    cur_ids[1] = 0;
    cur_ids[2]++;
    if (cur_ids[2] < neighbor_set.size()) return true;

    // If we made it this far, we are done.  Cycle neighbor back and return false.
    cur_ids[2] = 0;
    return false;
  }
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
    if (neighbors == 0 || neighbors > 8) {
      neighbors = SIZE;
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
  
  static size_t TestMulticell(emp::Random & random, const Config & config) {
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
          if (++orgs[pos] == config.threshold) {
            orgs[pos] = 1;

            size_t next_pos = RandomNeighbor(random, pos, config.neighbors);
            size_t & next_org = orgs[next_pos];

            // If the target is empty, put a new organism there.
            if (!next_org) {
              next_org = 1;
              num_orgs++;
            }

            // Otherwise if we don't restrain, reset organism there.
            else if (!config.restrain) {
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
                  const Config & config,
                  std::ostream & os=std::cout,
                  bool verbose=false,
                  bool headers=true)
  {
    World<CUR_SIZE,CUR_SIZE> world;
    
    // Only the first time through should we print the column headers.
    if (headers) {
      os << "world_x, world_y, " << config.GetHeaders();
      if (verbose) {
        for (size_t i=0; i < config.num_runs; i++) os << ", run" << i;
      }
      os << ", ave_time" << std::endl;
    }

    // Output current config info.
    os << CUR_SIZE << ", " << CUR_SIZE << ", " << config.AsCSV();
    
    double total_time = 0.0;
    for (size_t i = 0; i < config.num_runs; i++) {
      size_t cur_time = world.TestMulticell(random, config);
      if (verbose) os << ", " << cur_time;
      total_time += (double) cur_time;
    }
    os << ", " << (total_time / (double) config.num_runs) << std::endl;
  }

  // Run all of the configurations in an entire set.
  static void Run(emp::Random & random,
                  ConfigSet config_set,
                  std::ostream & os=std::cout,
                  bool verbose=false,
                  bool headers=true)
  {
    // Build a world of the correct size.
    if (CUR_SIZE <= config_set.max_size) {
      const size_t num_configs = config_set.GetSize();
      for (size_t i = 0; i < num_configs; i++) {
        Run(random, config_set.GetConfig(), std::cout, verbose, i==0 && headers);
        config_set.Next();
      }
      headers = false; // Don't do headers more than once.
    }

    // Do the recursive call.
    WorldSet<WORLD_SIZES...>::Run(random, config_set, os, verbose, headers);
  }

};

void PrintHelp(const std::string & name) {
  std::cout << "Format: " << name << " [OPTIONS...]\n"
            << "Options include:\n"
            << " -h         : This message (--help).\n"
            << " -n [SIZES] : Comma separated neighborhood sizes (--neighbors).\n"
            << "\nExample:  " << name << " -n 0,4,8\n"
            << std::endl;
}

int main(int argc, char* argv[])
{
  emp::vector<std::string> args = emp::cl::args_to_strings(argc, argv);
  emp::vector<size_t>  neighbor_set = {8};  // Default to size 8 neighborhoods.

  if (emp::Has<std::string>(args, "-h") || emp::Has<std::string>(args, "--help")) {
    PrintHelp(args[0]); exit(0);
  }

  size_t arg_id = 1;
  while (arg_id < args.size()) {
    const std::string cur_arg = args[arg_id++];
    if (cur_arg == "-n" || cur_arg == "--neighbors") {
      if (arg_id >= args.size()) { std::cout << "ERROR: Must provide neighborhood sizes!\n"; exit(1); }
      neighbor_set = emp::from_strings<size_t>(emp::slice(args[arg_id++], ','));
    }
  }

  emp::Random random;
  ConfigSet config_set;
  emp::Append(config_set.restrain_set, true, false);
  emp::Append(config_set.threshold_set, 4, 8, 16, 32);
  // emp::Append(config_set.neighbor_set, 0, 4, 6, 8);
  config_set.neighbor_set = neighbor_set;
  config_set.num_runs = 100;

  WorldSet<2,4,8,16,32,64,128>::Run(random, config_set, std::cout, false);
}
