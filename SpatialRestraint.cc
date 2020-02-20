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
  size_t cells_side = 32;    // How many cells are on a side of the (square) multi-cell?
  bool restrain = false;     // Should cells refrain from overwriting each other?
  size_t threshold = 16;     // How many resources are needed to produce an offspring?
  size_t neighbors = 8;      // Num neighbors to consider for offspring (0=well mixed; 4,6,8 => 2D)

  size_t num_runs = 100;     // How many runs should we do with the above configuration?
  bool verbose = false;      // Should we print timings for each replicate?

  size_t GetWidth() const { return cells_side; }
  size_t GetHeight() const { return cells_side; }
  size_t GetSize() const { return cells_side * cells_side; }

  std::string AsCSV() const {
    return emp::to_string(GetWidth(), ", ", GetHeight(), ", ", threshold, ", ", restrain, ", ", neighbors);
  }

  size_t ToPos(size_t x, size_t y) const { return x + y * cells_side; }
  size_t ToX(size_t pos) const { return pos % cells_side; }
  size_t ToY(size_t pos) const { return pos / cells_side; }
};

// Set of parameters to use.
struct ConfigSet {
  emp::vector<size_t> side_set{16};         ///< Values for cell side (16x16 = size 256 multicell.)
  emp::vector<bool> restrain_set{0,1};      ///< Values to use for restrain.
  emp::vector<size_t> threshold_set{16};    ///< Values to use for threhsold.
  emp::vector<size_t> neighbor_set{8};      ///< How many neighbors should we look at; 

  emp::array<size_t, 4> cur_ids{0,0,0,0};   ///< Which settings are we going to use next?

  size_t num_runs = 100;                    ///< How many times should we run each configuration?
  bool verbose = false;                     ///< Should we print data for each replicate?

  std::string GetHeaders() const { return "width, height, threshold, restrain, neighbors"; }

  size_t GetSize() const {
    return side_set.size() * restrain_set.size() * threshold_set.size() * neighbor_set.size();
  }

  Config GetConfig() const {
    return Config{ side_set[cur_ids[0]],
                   restrain_set[cur_ids[1]],
                   threshold_set[cur_ids[2]],
                   neighbor_set[cur_ids[3]],
                   num_runs,
                   verbose
                  };
  }

  bool Next() {
    // First try to toggle sides.
    cur_ids[0]++;
    if (cur_ids[0] < side_set.size()) return true;

    // Next, cycle sides back to the beginning and try to move to the next restrain.
    cur_ids[0] = 0;
    cur_ids[1]++;
    if (cur_ids[1] < restrain_set.size()) return true;

    // Next, cycle restrain back to the beginning and try to move to the next threshold.
    cur_ids[1] = 0;
    cur_ids[2]++;
    if (cur_ids[2] < threshold_set.size()) return true;

    // If we're done with thresholds, cycle back and try neighbors.
    cur_ids[2] = 0;
    cur_ids[3]++;
    if (cur_ids[3] < neighbor_set.size()) return true;

    // If we made it this far, we are done.  Cycle neighbor back and return false.
    cur_ids[3] = 0;
    return false;
  }
};

// World to test a single multicell for how long it takes to replicate.

struct World {

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

  static size_t RandomNeighbor(emp::Random & random, size_t pos, const Config & config)
  {
    if (config.neighbors == 0 || config.neighbors > 8) {
      return random.GetUInt(config.GetSize());
    }

    const size_t x = config.ToX(pos);
    const size_t y = config.ToY(pos);
    size_t next_x = (size_t) -1;
    size_t next_y = (size_t) -1;

    while (next_x >= config.GetWidth() || next_y >= config.GetHeight()) {
      const size_t dir = random.GetUInt(config.neighbors);  // Direction for offspring.
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

    emp_assert(config.ToPos(next_x, next_y) < config.GetSize(), next_x, next_y, config.cells_side);

    return config.ToPos(next_x, next_y);
  }

  static void Print(const emp::vector<size_t> & mc, const Config & config) {
    emp_assert(mv.size() == config.GetSize());
    size_t pos = 0;
    for (size_t y = 0; y < config.GetHeight(); y++) {
      for (size_t x = 0; x < config.GetWidth(); x++) {
      	std::cout << " " << ToChar(mc[pos++]);
      }
      std::cout << std::endl;
    }
  }
  
  static size_t TestMulticell(emp::Random & random, const Config & config) {
    // Setup initial multicell to be empty; keep count of resources in each cell.
    const size_t mc_size = config.GetSize();
    emp::vector<size_t> orgs(mc_size, 0);
    size_t time = 0;
    
    // Inject a cell in the middle.
    const size_t start_pos = config.ToPos(config.GetWidth()/2, config.GetHeight()/2);
    orgs[start_pos] = 1; // Live cells start at one resource.
    size_t num_orgs = 1;    

    // Loop through updates until cell is full.
    while (num_orgs < mc_size) {
      // std::cout << "Update: " << time << "; num_orgs: " << num_orgs << std::endl;
      // if (time % 100 == 0) Print(config, orgs);

      time++;

      // Loop through cells, probabilistically incrementing resources.
      for (size_t pos = 0; pos < mc_size; pos++) {
        // See if this cell gets a unit of resource.
        if (orgs[pos] && random.P(0.5)) {
          // If so, give resource and see if it replicates.
          if (++orgs[pos] == config.threshold) {
            orgs[pos] = 1;

            size_t next_pos = RandomNeighbor(random, pos, config);
            size_t & next_org = orgs[next_pos];

            // If the target is empty, put a new organism there.
            if (!next_org) {
              next_org = 1;
              num_orgs++;
            }

            // Otherwise if we don't restrain, reset existing organism there.
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

void Run(emp::Random & random,
                const Config & config,
                std::ostream & os=std::cout)
{
  World world;
  
  // Output current config info.
  os << config.AsCSV();
  
  double total_time = 0.0;
  for (size_t i = 0; i < config.num_runs; i++) {
    size_t cur_time = world.TestMulticell(random, config);
    if (config.verbose) os << ", " << cur_time;
    total_time += (double) cur_time;
  }
  os << ", " << (total_time / (double) config.num_runs) << std::endl;
}

// Run all of the configurations in an entire set.
void Run(emp::Random & random,
                ConfigSet config_set,
                std::ostream & os=std::cout)
{
  // Start with column headers.
  os << config_set.GetHeaders();
  if (config_set.verbose) {
    for (size_t i=0; i < config_set.num_runs; i++) os << ", run" << i;
  }
  os << ", ave_time" << std::endl;

  // Build a world of the correct size.
  const size_t num_configs = config_set.GetSize();
  for (size_t i = 0; i < num_configs; i++) {
    Run(random, config_set.GetConfig(), std::cout);
    config_set.Next();
  }
}

void PrintHelp(const std::string & name) {
  std::cout << "Format: " << name << " [OPTIONS...]\n"
            << "Options include:\n"
            << " -c [SIDE_SIZES] : Cells on side of (square) multicell (--cells_side) [16]\n"
            << " -d [COUNT]      : How many data replicates should we run? (--data_count) [100]\n"
            << " -h              : This message (--help).\n"
            << " -n [SIZES]      : Comma separated neighborhood sizes (--neighbors) [8].\n"
            << " -r [RESTRAINS]  : Should cells restrain? (--restrains) [0,1].\n"
            << " -t [THRESHOLDS] : Comma separated cell-repro thresholds (--thresholds) [16].\n"
            << " -v              : Use verbose data printing ALL results (--verbose) [false]"
            << "\nExample:  " << name << " -n 0,4,8 -r 0,1 -t 4,8,16,32 -d 100\n"
            << std::endl;
}

void ProcessCommandLine(ConfigSet & config_set, int argc, char* argv[]) {
  emp::vector<std::string> args = emp::cl::args_to_strings(argc, argv);

  if (emp::Has<std::string>(args, "-h") || emp::Has<std::string>(args, "--help")) {
    PrintHelp(args[0]); exit(0);
  }

  size_t arg_id = 1;
  while (arg_id < args.size()) {
    const std::string cur_arg = args[arg_id++];

    if (cur_arg == "-c" || cur_arg == "--cells_side") {
      if (arg_id >= args.size()) { std::cout << "ERROR: Must provide side-lengths to use!\n"; exit(1); }
      config_set.side_set = emp::from_strings<size_t>(emp::slice(args[arg_id++], ','));
    }

    else if (cur_arg == "-d" || cur_arg == "--data_count") {
      if (arg_id >= args.size()) { std::cout << "ERROR: Must provide data count!\n"; exit(1); }
      config_set.num_runs = emp::from_string<size_t>(args[arg_id++]);
    }

    else if (cur_arg == "-n" || cur_arg == "--neighbors") {
      if (arg_id >= args.size()) { std::cout << "ERROR: Must provide neighborhood sizes!\n"; exit(1); }
      config_set.neighbor_set = emp::from_strings<size_t>(emp::slice(args[arg_id++], ','));
    }

    else if (cur_arg == "-r" || cur_arg == "--restrains") {
      if (arg_id >= args.size()) { std::cout << "ERROR: Must provide restrain values!\n"; exit(1); }
      config_set.restrain_set = emp::from_strings<bool>(emp::slice(args[arg_id++], ','));
    }

    else if (cur_arg == "-t" || cur_arg == "--thresholds") {
      if (arg_id >= args.size()) { std::cout << "ERROR: Must provide threshold values!\n"; exit(1); }
      config_set.threshold_set = emp::from_strings<size_t>(emp::slice(args[arg_id++], ','));
    }

    else if (cur_arg == "-v" || cur_arg == "--verbose") {
      config_set.verbose = true;
    }

    else {
      std::cerr << "ERROR: Unknown option " << cur_arg << "\n";
      exit(1);
    }
  }
}

int main(int argc, char* argv[])
{
  emp::Random random;
  ConfigSet config_set;
  ProcessCommandLine(config_set, argc, argv);
  Run(random, config_set, std::cout);
}
