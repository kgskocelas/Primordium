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

#include "config/SettingConfig.h"
#include "tools/BitVector.h"
#include "tools/Distribution.h"
#include "tools/StreamManager.h"
#include "tools/string_utils.h"
#include "tools/vector_utils.h"

#include "Multicell.h"

/// Information about a full multi-cell organism
struct Organism {
  size_t num_ones = 0;
  double gen = 0.0;
  double repro_time = 0.0;

  Organism(size_t in_ones, double in_gen=0.0, double in_repro_time=0.0)
   : num_ones(in_ones), gen(in_gen), repro_time(in_repro_time) { }
  Organism(const Organism &) = default;
  Organism & operator=(const Organism &) = default;
};

struct Population {
  emp::vector<Organism> orgs;        ///< Actual organisms in this population.
  size_t num_samples;                ///< Number of samples used to approximate repro distributions.
  emp::TimeQueue<size_t> org_queue;  ///< Track times for when orgs will replicate.
  double ave_gen = 0.0;              ///< Current generation of population (ave across orgs)

  /// We need to store the time distribution for reproduction.
  /// The outer vector is the number of ones.  The inner vector is a set of
  /// how long multicells took to replicate with that number of ones.
  /// Dimensions are GENOME_SIZE+1 -by- NUM_SAMPLES
  emp::vector< emp::vector<double> > repro_cache;  

  // Shared resources with Experiment
  Multicell & multicell;
  emp::Random & random;
  emp::StreamManager & stream_manager;

  Population(size_t pop_size, size_t ancestor_1s, size_t _samples,
             Multicell & _mc, emp::Random & _rand, emp::StreamManager & _smanager)
    : orgs(pop_size, ancestor_1s), num_samples(_samples)
    , repro_cache(_mc.genome_size + 1)
    , multicell(_mc), random(_rand), stream_manager(_smanager)
  {
  }

  void Reset(size_t pop_size, size_t ancestor_1s, bool reset_cache=true) {
    orgs.resize(0, ancestor_1s);
    orgs.resize(pop_size, ancestor_1s);
    org_queue.Reset();
    ave_gen = 0;
    if (reset_cache) {
      repro_cache.resize(0);
      repro_cache.resize(multicell.genome_size + 1);
    }
  }

  double CalcAveOnes() {
    double total_bits = 0.0;
    for (Organism & org : orgs) total_bits += (double) org.num_ones;
    return total_bits / (double) orgs.size();
  }

  double CalcAveGen() {
    double total_gen = 0.0;
    for (Organism & org : orgs) total_gen += org.gen;
    return total_gen / (double) orgs.size();
  }

  Organism CalcAveOrg() {
    Organism total_org(0);
    for (Organism & org : orgs) {
      total_org.num_ones += (double) org.num_ones;
      total_org.gen += org.gen;
      total_org.repro_time += org.repro_time;
    }
    return Organism(total_org.num_ones/orgs.size(),
                    total_org.gen / (double) orgs.size(),
                    total_org.repro_time / (double) orgs.size());
  }

  double CalcReproDuration(size_t num_ones) {
    emp::vector<double> & cur_cache = repro_cache[num_ones];
    size_t sample_id = random.GetUInt(num_samples);
    if (sample_id < cur_cache.size()) return cur_cache[sample_id];

    multicell.start_1s = num_ones;
    multicell.SetupConfig();
    multicell.InjectCell(multicell.MiddlePos());
    double run_time = multicell.Run().GetReproTime();
    // std::cout << "run_time = " << run_time << std::endl;

    cur_cache.push_back(run_time);
    return run_time;
  }

  double CalcBirthTime(size_t num_ones) {
    return CalcReproDuration(num_ones) + org_queue.GetTime();
  }

  double CalcAveReproDuration() {
    double total_rt = 0.0;
    for (Organism & org : orgs) total_rt += CalcReproDuration(org.num_ones);
    return total_rt / (double) orgs.size();
  }

  void NextBirth() {
    size_t parent_id = org_queue.Next();
    Organism & parent = orgs[parent_id];

    // std::cout << "DEBUG: NextBirth with parent_id=" << parent_id << std::endl;

    // If this organism has updated, skip it.
    if (parent.repro_time != org_queue.GetTime()) {
      // std::cout << "DEBUG: ...skipped; parent.repro_time=" << parent.repro_time
      //           << ", but real time=" << org_queue.GetTime() << std::endl;
      return;
    }

    // Figure out where the offspring would go.
    size_t offspring_id = random.GetUInt(orgs.size());
    Organism & offspring = orgs[offspring_id];

    // std::cout << "DEBUG: ...offspring_id=" << offspring_id << std::endl;

    ave_gen -= offspring.gen / (double) orgs.size();      // Remove old org from gen average.
    if (parent_id != offspring_id) {                      // If the parent is not being replaced...
      offspring = parent;                                 //   copy parent to offspring.
      parent.repro_time = CalcBirthTime(parent.num_ones); //   figure out parent's NEXT repro time.
      org_queue.Insert(parent_id, parent.repro_time);     //   schedule parent for next repro
    }
    offspring.gen += 1.0;                                 // Update offspring's generation.
    ave_gen += offspring.gen / (double) orgs.size();      // Add new org to gen average.

    // Handle mutations in the offspring.
    if (random.P(multicell.mut_prob)) {
      double prob1 = ((double) offspring.num_ones) / (double) multicell.genome_size;
      if (random.P(prob1)) offspring.num_ones--;
      else offspring.num_ones++;
    }


    // Schedule offspring to give birth.
    offspring.repro_time = CalcBirthTime(offspring.num_ones);
    org_queue.Insert(offspring_id, offspring.repro_time);
  }

  void Run(double max_gen, const std::string run_name="", bool verbose=false) {
    // Setup the time queue.
    for (size_t i = 0; i < orgs.size(); i++) {
      double repro_time = CalcBirthTime(orgs[i].num_ones);
      org_queue.Insert(i, repro_time);
      orgs[i].repro_time = repro_time;
    }

    // If verbose or print_reps is turned on, we need to track the current generation.
    if (verbose || run_name.size()) {
      std::ostream & os(stream_manager.get_ostream(run_name));

      bool print_both = verbose && run_name.size();  // Should we send output to both places?

      os << "#generation, ave_ones, ave_repro_time\n";
      if (print_both) std::cout << "#generation, ave_ones, ave_repro_time\n";

      double next_gen = -1.0;
      std::string out_line = "";
      while (ave_gen < max_gen) {
        if (ave_gen > next_gen) {
          next_gen += 1.0;
          out_line = emp::to_string((size_t) next_gen,
                                    ", ", CalcAveOnes(),
                                    ", ", CalcAveReproDuration()
                                   );
          os << out_line << std::endl;
          if (print_both) std::cout << out_line << std::endl;
        }
        NextBirth();
      }
    }

    else {
      while (ave_gen < max_gen) {
        NextBirth();
      }
    }
  }

  void PrintData(std::ostream & os=std::cout) {
    // Count up the number of organism with each bit count.
    emp::vector<size_t> bit_counts(repro_cache.size(), 0);
    for (Organism & org : orgs) bit_counts[org.num_ones]++;

    // And print the results.
    for (size_t i = 0; i < bit_counts.size(); i++) {
      os << ", " << bit_counts[i];
    }
    os << ", " << CalcAveOnes();
  }
};

struct Experiment {
  emp::Random random;
  emp::SettingConfig config;
  std::string exe_name;      ///< Name of executable used to start this run.
  Multicell multicell;

  size_t gen_count = 0;          ///< Num generations to evolve (zero for analyze multicells)
  size_t pop_size = 200;         ///< Num organisms in the population.
  size_t sample_size = 100;      ///< Num multicells to sample for each genotype.
  bool balance_predict = false;  ///< Try to predict the mutation-selection balance.
  bool print_reps = false;       ///< Should we print results for every replicate?
  bool print_trace = false;      ///< Should we show each step of a multicell?
  bool reset_cache = false;      ///< Share the cache by default.
  bool verbose = false;          ///< Should we print extra information during the run?

  emp::StreamManager stream_manager;  ///< Manage files
  std::string evolution_filename;     ///< Output filename for evolution summary data.
  std::string multicell_filename;     ///< Output filename for multicell summary data.

  using TreatmentResults = emp::vector<RunResults>;
  using MulticellResults = emp::vector<TreatmentResults>;

  MulticellResults base_results;

  Experiment(emp::vector<std::string> & args) : multicell(random) {
    exe_name = args[0];

    // Setup all command-line options that the system should use.  In general, lower-case
    // letters are used to control model parameters, while capital letters are used to control
    // output.  The one exception is -h for '--help' which is otherwise too standard.
    // The order below sets the order that combinations are tested in. 
    // AVAILABLE OPTION FLAGS: efjklqwxyz ACDFGHIJKLNOQRSUVWXYZ

    config.AddComboSetting<size_t>("data_count", "Number of times to replicate each run", 'd') = { 100 };
    config.AddComboSetting("ancestor_1s", "How many 1s in starting cell?", 'a',
                           multicell.start_1s, "NumOnes...") = { 50 };
    config.AddComboSetting("unrestrained_cost", "Per-cell cost for unrestrained", 'u',
                           multicell.unrestrained_cost, "Costs...") = { 0.0 };
    config.AddComboSetting("mut_prob",   "Probability of mutation in offspring", 'm',
                           multicell.mut_prob, "Probs...") = { 0.0 };
    config.AddComboSetting("time_range", "Rep time = 100.0 + random(time_range)", 't',
                           multicell.time_range, "TimeUnits...") = { 50 };
    config.AddComboSetting("neighbors",  "Neighborhood size for replication", 'n',
                           multicell.neighbors, "Sizes...") = { 8 };
    config.AddComboSetting("restrain",   "Num ones in genome for restraint?", 'r',
                           multicell.restrain, "NumOnes...") = { 50 };
    config.AddComboSetting("bit_size",   "Number of bits in genome?", 'b',
                           multicell.genome_size, "NumBits...") = { 100 };
    config.AddComboSetting("cells_side", "Cells on side of (square) multicell", 'c',
                           multicell.cells_side, "NumCells...") = { 32 };

    config.AddAction("one_check", "Make restrained check only one cell to find empty.", 'o',
                     [this](){ multicell.one_check = true; } );

    config.AddSetting("gen_count",   "Num generations to evolve (0=analyze only)", 'g',
                      gen_count, "NumGens") = { 0 };
    config.AddSetting("pop_size",    "Number of organisms in the population.", 'p',
                      pop_size, "NumOrgs") = { 200 };
    config.AddSetting("sample_size", "Num multicells sampled for distributions.", 's',
                      sample_size, "NumSamples") = { 200 };
                      

    config.AddAction("balance_predict", "Predict the mutation-selection balance [NOT YET IMPLEMENTED!]", 'B',
                     [this](){ balance_predict = true; } );
    config.AddAction("help", "Print full list of options", 'h',
                     [this](){
                       config.PrintHelp(exe_name, " -n 0,4,8 -r 0,1 -t 4,8,16,32 -d 100");
                       exit(1);
                      } );
    config.AddSetting("evolution_filename", "Filename for multicell data", 'E',
                      evolution_filename, "Filename") = "evolution.dat";
    config.AddAction("independent_caches", "Use a distinct cache for each run", 'i',
                     [this](){ reset_cache = true; } );
    config.AddSetting("multicell_filename", "Filename for multicell data", 'M',
                      multicell_filename, "Filename") = "multicell.dat";
    config.AddAction("print_reps", "Print data for each replicate", 'P',
                     [this](){ print_reps = true; } );
    config.AddAction("trace", "Show each step of replicates (multicell or population)", 'T',
                     [this](){ print_trace = true; } );
    config.AddAction("verbose", "Print extra information during the run", 'v',
                     [this](){ verbose = true; } );

    // Process the command-line options
    config.ProcessOptions(args);

    // Fail if there are any unknown args.
    if (config.HasUnusedArgs()) {
      std::cerr << "ERROR: Unknown options: " << emp::to_string(config.GetUnusedArgs()) << "\n";
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
    const size_t num_runs = config.GetValue<size_t>("data_count");
    const size_t combo_id = config.GetComboID();
    TreatmentResults & treatment_results = base_results[combo_id];
    treatment_results.resize(num_runs);

    // Conduct all replicates and output the information.    
    for (size_t i = 0; i < num_runs; i++) {
      treatment_results[i] = TestMulticell();
      if (print_reps) os << ", " << treatment_results[i].GetReproTime();
    }

    return treatment_results;
  }

  RunResults SummarizeTreatment(std::ostream & os=std::cout) {
    const size_t num_runs = config.GetValue<size_t>("data_count");
    const size_t combo_id = config.GetComboID();

    // Setup room for the data being collected.
    TreatmentResults & treatment_results = base_results[combo_id];
    treatment_results.resize(num_runs);

    // Conduct all replicates and output the information.    
    RunResults total_results(multicell.genome_size);
    for (size_t i = 0; i < num_runs; i++) {
      if (verbose) std::cout << " ... run " << i << std::endl;
      treatment_results[i] = TestMulticell();
      if (print_reps) os << ", " << treatment_results[i].GetReproTime();
      total_results += treatment_results[i];
    }

    return total_results /= (double) num_runs;
  }

  /// Given the current configuration options, evolve a set of runs.
  void EvolveTreatment(std::ostream & os=std::cout) {
    const size_t num_runs = config.GetValue<size_t>("data_count");
    const size_t num_samples = config.GetValue<size_t>("sample_size");
    const size_t pop_size = config.GetValue<size_t>("pop_size");
    const size_t ancestor_1s = config.GetValue<size_t>("ancestor_1s");
    const size_t gen_count = config.ComboValues<size_t>("gen_count")[0];

    Population pop(pop_size, ancestor_1s, num_samples, multicell, random, stream_manager);
    for (size_t run_id = 0; run_id < num_runs; run_id++) {
      std::cout << "START Treatment #" << config.GetComboID()
                << " : Run " << run_id << std::endl;
      std::string run_name =
        print_trace ? emp::to_string('t',config.GetComboID(),'r',run_id,".dat") : "";
      pop.Reset(pop_size, ancestor_1s, reset_cache);
      pop.Run(gen_count, run_name, verbose);

      os << config.CurComboString(", ");  // Output current setting combination data.
      pop.PrintData(os);             // Output data for THIS population.
      os << std::endl;
    }
  }

  /// Step through all configurations and collect multicell data for each.
  void RunMulticells(std::ostream & os) {
    // Print column headers.
    os << "#" << config.GetComboHeaders();
    if (print_reps) {
      const size_t num_runs = config.GetValue<size_t>("data_count");
      for (size_t i=0; i < num_runs; i++) os << ", run" << i;
    }
    os << ", ave_time, frac_restrain" << std::endl;

    // Setup the correct collection for the treatments.
    base_results.resize(config.CountCombos());

    // Loop through configuration combonations to test.
    config.ResetCombos();
    do {
      std::cout << "START Treatment #" << config.GetComboID()
                << " / " << base_results.size()
                << std::endl
                << "  " << config.CurComboString(", ", true, true)
                << std::endl;

      os << config.CurComboString(", ");  // Output current setting combination data.

      RunResults treatment_results = SummarizeTreatment(os);

      os << ", " << treatment_results.GetReproTime()
         << ", " << (treatment_results.CountRestrained(multicell.restrain) / (double) multicell.GetSize())
         << std::endl;
    } while (config.NextCombo());
  }

  void RunEvolution(std::ostream & os) {
    // Print column headers.
    os << "#" << config.GetComboHeaders();
    const size_t max_bits = config.MaxValue<size_t>("bit_size");
    for (size_t i=0; i <= max_bits; i++) os << ", " << i << "-ones";
    os << ", ave_ones" << std::endl;

    config.ResetCombos();
    do {
      EvolveTreatment(os);
    } while (config.NextCombo());
  }

  // Run all of the configurations in an entire set.
  void Run() {
    size_t gen_count = config.GetValue<size_t>("gen_count");
    std::string evolution_filename = config.GetValue<std::string>("evolution_filename");
    std::string multicell_filename = config.GetValue<std::string>("multicell_filename");

    // If we have a generation count, collect evolution data.
    if (gen_count) RunEvolution(stream_manager.get_ostream(evolution_filename));

    // Otherwise collect information on multicells.
    else RunMulticells(stream_manager.get_ostream(multicell_filename));
  }
};

#endif