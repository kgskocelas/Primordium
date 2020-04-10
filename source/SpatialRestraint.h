/**
 *  @copyright Copyright (C) Michigan State University, MIT Software license; see doc/LICENSE.md
 *  @date 2020.
 *
 *  @file  SpatialRestraint.h
 *  @brief Manages analysis of multicells and how they handle cells that exhibit restraint.
 *  @note Status: BETA
 *
 *  Cells are either RESTRAINED or UNRESTRAINED
 *  Organisms are "multi-cells" that start as a single cell and automatically replicate when full.
 *  Cells replicate within the multicell when they have collected enough resource.
 *  - Restrained cells replicate into an empty neighbor position (or fail if no empty)
 *  - Unrestrained cells replicate into a random position, regardless of if anything is there.
 */

#ifndef SPATIAL_RESTRAINT_H
#define SPATIAL_RESTRAINT_H

#include <iostream>
#include <fstream>
#include <set>

#include "config/command_line.h"
#include "config/SettingCombos.h"
#include "tools/BitVector.h"
#include "tools/Distribution.h"
#include "tools/StreamManager.h"
#include "tools/string_utils.h"
#include "tools/vector_utils.h"

#include "Multicell.h"

/// Information about a full multi-cell organism
struct Organism {
  size_t num_ones = 0;
  size_t gen = 0;

  Organism(size_t in_ones) : num_ones(in_ones) { }
  Organism(const Organism &) = default;
  Organism & operator=(const Organism &) = default;
};

struct Population {
  emp::vector<Organism> orgs;
  size_t num_samples;
  Multicell & multicell;
  emp::vector< emp::vector<double> > repro_cache;

  Population(size_t pop_size, size_t initial_1s, size_t _samples, Multicell & _mc)
    : orgs(pop_size, initial_1s), num_samples(_samples), multicell(_mc)
    , repro_cache(multicell.genome_size + 1)
  {
    for (auto & size_cache : repro_cache) size_cache.resize(num_samples, 0.0);
  }
};

struct Experiment {
  emp::Random random;
  emp::SettingCombos combos;
  std::string exe_name;      ///< Name of executable used to start this run.
  Multicell multicell;

  size_t gen_count = 0;      ///< Num generations to evolve (zero for analyze multicells)
  size_t pop_size = 200;     ///< Num organisms in the population.
  size_t sample_size = 100;  ///< Num multicells to sample for each genotype.
  bool print_reps = false;   ///< Should we print results for every replicate?
  bool print_trace = false;  ///< Should we show each step of a multicell?

  std::string evolution_filename;  ///< Output filename for evolution summary data.
  std::string multicell_filename;  ///< Output filename for multicell summary data.

  using TreatmentResults = emp::vector<RunResults>;
  using MulticellResults = emp::vector<TreatmentResults>;

  MulticellResults base_results;

  Experiment(emp::vector<std::string> & args) : multicell(random) {
    exe_name = args[0];

    combos.AddSetting("time_range", "Rep time = 100.0 + random(time_range)", 't',
                       multicell.time_range, "TimeUnits...") = { 50 };
    combos.AddSetting("neighbors",  "Neighborhood size for replication", 'n',
                      multicell.neighbors, "Sizes...") = { 8 };
    combos.AddSetting("cells_side", "Cells on side of (square) multicell", 'c',
                      multicell.cells_side, "NumCells...") = { 32 };
    combos.AddSetting("bit_size",   "Number of bits in genome?", 'b',
                      multicell.genome_size, "NumBits...") = { 10 };
    combos.AddSetting("restrain",   "Num ones in genome for restraint?", 'r',
                      multicell.restrain, "NumOnes...") = { 5 };
    combos.AddSetting("initial_1s", "How many 1s in starting cell?", 'i',
                      multicell.start_1s, "NumOnes...") = { 5 };
    combos.AddSetting("mut_prob",   "Probability of mutation in offspring", 'm',
                      multicell.mut_prob, "Probs...") = { 0.0 };
    combos.AddSetting<size_t>("data_count", "Number of times to replicate each run", 'd') = { 100 };

    combos.AddSingleSetting("gen_count",   "Num generations to evolve (0=analyze only)", 'g',
                      gen_count, "NumGens") = { 0 };
    combos.AddSingleSetting("pop_size",    "Number of organisms in the population.", 'p',
                      pop_size, "NumOrgs") = { 200 };
    combos.AddSingleSetting("sample_size", "Num. multicells sampled for fitness distribution.", 's',
                      sample_size, "NumSamples") = { 200 };
                      

    combos.AddAction("help", "Print full list of options", 'h',
                     [this](){
                       combos.PrintHelp(exe_name, " -n 0,4,8 -r 0,1 -t 4,8,16,32 -d 100");
                       exit(1);
                      } );
    combos.AddAction("print_reps", "Print timings for each replicate", 'P',
                     [this](){ print_reps = true; } );
    combos.AddAction("trace", "Show each step of a multicell", 'T',
                     [this](){ print_trace = true; } );
    combos.AddSingleSetting("evolution_filename", "Filename for multicell data", 'E',
                            evolution_filename, "Filename") = { "evolution.dat" };
    combos.AddSingleSetting("multicell_filename", "Filename for multicell data", 'M',
                            multicell_filename, "Filename") = { "multicell.dat" };

    // Process the command-line options
    args = combos.ProcessOptions(args);

    // Fail if there are any unknown args.
    if (args.size()) {
      std::cerr << "ERROR: Unknown options: " << emp::to_string(args) << "\n";
      exit(2);
    }
  }


  RunResults TestMulticell() {
    multicell.SetupConfig();

    // Inject a cell in the middle.
    const size_t start_pos = multicell.MiddlePos();
    multicell.InjectCell(start_pos);

    // Do the run!
    return multicell.Run(print_trace);
  }

  TreatmentResults & RunTreatment(std::ostream & os=std::cout) {
    const size_t num_runs = combos.GetValue<size_t>("data_count");
    const size_t combo_id = combos.GetComboID();
    TreatmentResults & treatment_results = base_results[combo_id];
    treatment_results.resize(num_runs);

    // Conduct all replicates and output the information.    
    for (size_t i = 0; i < num_runs; i++) {
      treatment_results[i] = TestMulticell();
      if (print_reps) os << ", " << treatment_results[i].run_time;
    }

    return treatment_results;
  }

  RunResults SummarizeTreatment(std::ostream & os=std::cout) {
    const size_t num_runs = combos.GetValue<size_t>("data_count");
    const size_t combo_id = combos.GetComboID();

    // Setup room for the data being collected.
    TreatmentResults & treatment_results = base_results[combo_id];
    treatment_results.resize(num_runs);

    // Conduct all replicates and output the information.    
    RunResults total_results(multicell.genome_size);
    for (size_t i = 0; i < num_runs; i++) {
      treatment_results[i] = TestMulticell();
      if (print_reps) os << ", " << treatment_results[i].run_time;
      total_results += treatment_results[i];
    }

    return total_results /= (double) num_runs;
  }

  /// Given the current configuration options, evolve a set of runs.
  void EvolveTreatment(std::ostream & os=std::cout) {
    // Setup the results cache.
    const size_t num_samples = combos.GetValue<size_t>("sample_size");
    const size_t pop_size = combos.GetValue<size_t>("pop_size");
    const size_t initial_1s = combos.GetValue<size_t>("initial_1s");

    Population pop(pop_size, initial_1s, num_samples, multicell);

    // Run through the generations.
    for (size_t gen=0; gen <= gen_count; gen++) {

    }
  }

  /// Step through all configurations and collect multicell data for each.
  void RunMulticells(std::ostream & os) {
    // Print column headers.
    os << combos.GetHeaders();
    if (print_reps) {
      const size_t num_runs = combos.GetValue<size_t>("data_count");
      for (size_t i=0; i < num_runs; i++) os << ", run" << i;
    }
    os << ", ave_time, frac_restrain" << std::endl;

    // Setup the correct collection for the treatments.
    base_results.resize(combos.CountCombos());

    // Loop through configuration combonations to test.
    combos.Reset();
    do {
      os << combos.CurString(", ");  // Output current setting combination data.

      RunResults treatment_results = SummarizeTreatment(os);

      os << ", " << treatment_results.run_time
         << ", " << (treatment_results.CountRestrained(multicell.restrain) / (double) multicell.GetSize())
         << std::endl;
    } while (combos.Next());
  }

  void RunEvolution(std::ostream & os) {
    combos.Reset();
    do {
      EvolveTreatment(os);
    } while (combos.Next());
  }

  // Run all of the configurations in an entire set.
  void Run() {
    size_t gen_count = combos.Values<size_t>("gen_count")[0];
    std::string evolution_filename = combos.Values<std::string>("evolution_filename")[0];
    std::string multicell_filename = combos.Values<std::string>("multicell_filename")[0];

    emp::StreamManager stream_manager;

    // If we have a generation count, collect evolution data.
    if (gen_count) RunEvolution(stream_manager.get_ostream(evolution_filename));

    // Otherwise collect information on multicells.
    else RunMulticells(stream_manager.get_ostream(multicell_filename));
  }
};

#endif