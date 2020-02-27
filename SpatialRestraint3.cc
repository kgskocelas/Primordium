#include <iostream>
#include <fstream>
#include "base/vector.h"
#include "config/command_line.h"
#include "config/SettingCombos.h"
#include "tools/Bool.h"
#include "tools/BitVector.h"
#include "tools/Distribution.h"
#include "tools/Random.h"
#include "tools/string_utils.h"
#include "tools/vector_utils.h"

/// Information about a single organism.
struct Organism {
  size_t resources = 0;  ///< How many resources have been collected?
  size_t threshold;      ///< How many resources are needed to replicate?
  emp::BitVector bits;   ///< What is underlying bit sequence?
};

/// Information needed to configure a run.
struct Config {
  size_t cells_side = 32;  ///< How many cells are on a side of the (square) multi-cell?
  size_t threshold = 16;   ///< How many resources are needed to produce an offspring?
  size_t neighbors = 8;    ///< Num neighbors to consider for offspring (0=well mixed; 4,6,8 => 2D)
  size_t bit_size = 10;    ///< How many bits in genome?
  size_t restrain = 5;     ///< How many ones in bit sequence for restraint?
  size_t start_1s = 5;     ///< How many ones in the starting organism?
  double mut_prob = 0.0;   ///< Probability of an offspring being mutated.

  Organism default_org;

  bool Set(emp::SettingCombos & combos) {
    cells_side = combos.GetValue<size_t>("cells_side");
    threshold = combos.GetValue<size_t>("threshold");
    neighbors = combos.GetValue<size_t>("neighbors");
    bit_size = combos.GetValue<size_t>("bit_size");
    restrain = combos.GetValue<size_t>("restrain");
    start_1s = combos.GetValue<size_t>("start_1s");
    mut_prob = combos.GetValue<double>("mut_prob");

    // If any setting combinations are impossible, reject.
    if (start_1s > bit_size) return false;

    // Setup the default organism.
    default_org.bits.Resize(bit_size);
    for (size_t i = 0; i < start_1s; i++) default_org.bits[i] = 1;

    return true;
  }

  size_t GetWidth() const { return cells_side; }
  size_t GetHeight() const { return cells_side; }
  size_t GetSize() const { return cells_side * cells_side; }

  size_t ToPos(size_t x, size_t y) const { return x + y * cells_side; }
  size_t ToX(size_t pos) const { return pos % cells_side; }
  size_t ToY(size_t pos) const { return pos / cells_side; }
};

struct World {
  emp::Random random;
  emp::SettingCombos combos;
  emp::vector<Organism> orgs;
  Config config;
  bool verbose = false;

  emp::NegativeBinomial threshold_dist;

  World(int argc, char* argv[]) : threshold_dist(0.5, 16) {
    combos.AddSetting<size_t>("threshold",  "Resources needed to replicate", "-t") = { 16 };
    combos.AddSetting<size_t>("neighbors",  "Neighborhood size for replication", "-n") = { 8 };
    combos.AddSetting<size_t>("cells_side", "Cells on side of (square) multicell", "-c") = { 16 };
    combos.AddSetting<size_t>("data_count", "Number of times to replicate each run", "-d") = { 100 };
    combos.AddSetting<size_t>("bit_size",   "Number of bits in genome?", "-b") = { 10 };
    combos.AddSetting<size_t>("restrain",   "Num ones in genome for restraint?", "-r") = { 5 };
    combos.AddSetting<size_t>("start_1s",   "How many 1s in starting organism?", "-1") = { 5 };
    combos.AddSetting<double>("mut_prob",   "Probability of mutation in offspring", "-m") = { 0.0 };

    emp::vector<std::string> args = emp::cl::args_to_strings(argc, argv);

    if (emp::Has<std::string>(args, "-h") || emp::Has<std::string>(args, "--help")) {
      PrintHelp(args[0]); exit(0);
    }

    args = combos.ProcessOptions(args);

    // Scan through remaining args...
    size_t arg_id = 1;
    while (arg_id < args.size()) {
      const std::string cur_arg = args[arg_id++];

      if (cur_arg == "-v" || cur_arg == "--verbose") { verbose = true; }
      else {
        std::cerr << "ERROR: Unknown option " << cur_arg << "\n";
        exit(1);
      }
    }
  }

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

  size_t RandomNeighbor(size_t pos)
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

  // Print current resource levels in population.
  void Print() {
    emp_assert(orgs.size() == config.GetSize());
    size_t pos = 0;
    for (size_t y = 0; y < config.GetHeight(); y++) {
      for (size_t x = 0; x < config.GetWidth(); x++) {
      	std::cout << " " << ToChar(orgs[pos++].resources);
      }
      std::cout << std::endl;
    }
  }

  void ResetParent(Organism & org) {
    org.resources = 1;                                 // Reset resources
    org.threshold = threshold_dist.PickRandom(random); // Determine new threshold.
    org.bits.Resize(config.bit_size);                  // Make sure its correct size.
  }

  void DoBirth(Organism & offspring, Organism & parent, bool do_mutations=true) {
    offspring.resources = 1;
    offspring.threshold = threshold_dist.PickRandom(random);
    offspring.bits = parent.bits;
    if (do_mutations && random.P(config.mut_prob)) {
      size_t pos = random.GetUInt(offspring.bits.size()); // Find bit to mutate.
      offspring.bits[pos] ^= 1;                           // Toggle it!
    }
  }
  
  size_t TestMulticell() {
    // Setup initial multicell to be empty; keep count of resources in each cell.
    const size_t mc_size = config.GetSize();
    orgs.resize(0);           // Clear out any current organisms.
    orgs.resize(mc_size);  // Put in new organsims initialized to 0.
    size_t time = 0;
    
    // Inject a cell in the middle.
    const size_t start_pos = config.ToPos(config.GetWidth()/2, config.GetHeight()/2);
    Organism & inject_org = orgs[start_pos]; // Find the org position to inject.
    inject_org = config.default_org;         // Initialize injection to proper default;
    ResetParent(inject_org);                 // Do any extra setup for this organism.
    size_t num_orgs = 1;                     // Count current organsims!

    // Loop through updates until cell is full.
    while (num_orgs < mc_size) {
      // std::cout << "Update: " << time << "; num_orgs: " << num_orgs << std::endl;
      // if (time % 100 == 0) Print();

      time++;

      // Loop through cells, probabilistically incrementing resources.
      for (size_t pos = 0; pos < mc_size; pos++) {
        Organism & org = orgs[pos];

        // See if this cell gets a unit of resource.
        if (org.resources) {
          // If so, give resource and see if it replicates.
          if (++org.resources == org.threshold) {
            Organism & parent = org;
            ResetParent(parent);

            size_t next_pos = RandomNeighbor(pos);
            Organism & next_org = orgs[next_pos];

            // If the target is empty, put a new organism there.
            if (!next_org.resources) {
              DoBirth(next_org, parent);
              num_orgs++;
            }

            // Otherwise if we don't restrain, reset existing organism there.
            else if (org.bits.CountOnes() < config.restrain) {
              DoBirth(next_org, parent);
            }
          }
        }
      }
    }

    return time;
  }

  // Run all of the configurations in an entire set.
  void Run(std::ostream & os=std::cout) {
    const size_t num_runs = combos.GetValue<size_t>("data_count");

    // Print column headers.
    os << combos.GetHeaders();
    if (verbose) {
      for (size_t i=0; i < num_runs; i++) os << ", run" << i;
    }
    os << ", ave_time" << std::endl;

    // Loop through configuration combonations to test.
    do {
      os << combos.CurString();  // Output current setting combination data.
      config.Set(combos);        // Store current setting combination in config object.      

      threshold_dist.Setup(0.5, config.threshold);

      // Conduct all replicates and output the information.    
      double total_time = 0.0;
      for (size_t i = 0; i < num_runs; i++) {
        size_t cur_time = TestMulticell();
        if (verbose) os << ", " << cur_time;
        total_time += (double) cur_time;
      }
      os << ", " << (total_time / (double) num_runs) << std::endl;
    } while (combos.Next());
  }

  void PrintHelp(const std::string & name) {
    std::cout << "Format: " << name << " [OPTIONS...]\n"
              << "Options include:\n"
              << " -b [NUM_BITS]   : How many bits should each organism have? (--bit_size) [10]\n"
              << " -c [SIDE_SIZES] : Cells on side of (square) multicell (--cells_side) [16]\n"
              << " -d [COUNT]      : How many data replicates should we run? (--data_count) [100]\n"
              << " -h              : This message (--help).\n"
              << " -m [MUT_RATE]   : Probability of mutation in offspring (--mut_prob) [0.0]"
              << " -n [SIZES]      : Comma separated neighborhood sizes (--neighbors) [8].\n"
              << " -r [RESTRAINT]  : How many ones to restrain? (--restrains) [5].\n"
              << " -t [THRESHOLDS] : Comma separated cell-repro thresholds (--thresholds) [16].\n"
              << " -v              : Use verbose data printing ALL results (--verbose) [false]\n"
              << " -1 [NUM_ONES]   : How many 1's in the starting organism? (--start_1s) [5]\n"
              << "\nExample:  " << name << " -n 0,4,8 -r 0,1 -t 4,8,16,32 -d 100\n"
              << std::endl;
  }
};

int main(int argc, char* argv[])
{
  World world(argc, argv);
  world.Run(std::cout);
}
