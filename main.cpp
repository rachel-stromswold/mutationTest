#include "main.h"

std::vector<unsigned> sample_bernoulli(unsigned n, double p, std::mt19937& generator, unsigned n_trials=1) {
  std::bernoulli_distribution bern(p);
  std::vector<unsigned> ret(n_trials, 0);
  for (unsigned i = 0; i < n_trials; ++i) {
    unsigned berVal = 0;
    for (unsigned j = 0; j < n; ++j) {
      if ( bern( generator ) ) {
        berVal = berVal | (0x01 << j);
      }
    }
    
    ret[i] = berVal;
  }
  return ret;
}

/*std::vector<unsigned> sample_binomial_shuffle(unsigned n, double p, std::mt19937& generator, unsigned n_trials=1) {
  //use the modified procedure, store each mutation mask in the binTest array
  std::binomial_distribution<unsigned> bin(n, p);
  std::vector<unsigned> ret(n_trials);
  for (unsigned i = 0; i < n_trials; ++i) {
    unsigned num_ones = bin(generator);
    unsigned val;

    for (unsigned j = 0; j < num_ones; ++j) {
      std::uniform_int_distribution<unsigned> unif(0, j);
      unsigned t = 1 << unif(generator);
      if ((val & t) != 0) {
        val = val | (1 << j);
      } else {
        val = val | t;
      }
    }
    ret[i] = val;
  }
  return ret;
}

std::vector<unsigned> sample_binomial(unsigned n, double p, std::mt19937& generator, unsigned n_trials=1) {
  unsigned invert = 0;
  if (p > 0.5) {
    invert=(2 << n) - 1;
    p=1-p;
  }
  //use the modified procedure, store each mutation mask in the binTest array
  std::binomial_distribution<unsigned> bin(n, p);
  std::vector<std::uniform_int_distribution<unsigned>>
  std::vector<unsigned> ret(n_trials);
  for (unsigned i = 0; i < n_trials; ++i) {
    unsigned num_ones = bin(generator);

    //minus 1 because we start indexing from 0
    std::uniform_int_distribution<unsigned> unif(0, choose(n, num_ones) - 1);
    unsigned binVal = get_bit_stream(n, num_ones, unif(generator));
    ret[i] = invert ^ binVal;
  }
  return ret;
}

std::vector<unsigned> sample_binomial_slow(unsigned n, double p, std::mt19937& generator, unsigned n_trials=1) {
  //use the modified procedure, store each mutation mask in the binTest array
  std::binomial_distribution<unsigned> bin(n, p);
  std::vector<unsigned> ret(n_trials);
  std::vector<unsigned> choosevals(n+1);
  for (unsigned i = 0; 2*i <= n; ++i) {
    choosevals[i] = choose(n, i);
    choosevals[n-i] = choosevals[i];
  }
  for (unsigned i = 0; i < n_trials; ++i) {
    unsigned num_ones = bin(generator);

    //minus 1 because we start indexing from 0
    std::uniform_int_distribution<unsigned> dist(0, choose(n, num_ones) - 1);
    unsigned binVal = get_bit_stream_slow(n, num_ones, dist(generator));
    ret[i] = binVal;
  }
  return ret;
}

std::vector<unsigned> sample_poisson(unsigned n, double p, std::mt19937& generator, unsigned n_trials=1) {
  double lambda = n*log(1-p);
  //use the modified procedure, store each mutation mask in the binTest array
  std::poisson_distribution<unsigned> poiss(lambda);
  std::uniform_int_distribution<unsigned> unif(0, n-1);
  std::vector<unsigned> ret(n_trials);
  for (unsigned i = 0; i < n_trials; ++i) {
    unsigned n_strings = poiss(generator);
    unsigned val = 0;
    for (unsigned j = 0; j < n_strings; ++j) {
      val = val | (1 << unif(generator));
    }
    ret[i] = val;
  }
  return ret;
}*/

struct OccurrenceCounter {
  unsigned bernoulli_occurrences = 0;
  unsigned poisson_occurrences = 0;
  unsigned binomial_old_occurrences = 0;
  unsigned binomial_new_occurrences = 0;
  unsigned digit_occurrences = 0;
};

struct TimingStats {
  unsigned bernoulli_total = 0;
  unsigned poisson_total = 0;
  unsigned binomial_old_total = 0;
  unsigned binomial_new_total = 0;
  double hybrid_poisson_total = 0;
  double hybrid_binomial_total = 0;

  double bernoulli_avg = 0;
  double poisson_avg = 0;
  double binomial_old_avg = 0;
  double binomial_new_avg = 0;
  double hybrid_poisson_avg = 0;
  double hybrid_binomial_avg = 0;
};

TimingStats test_non_hybrid(unsigned n_trials, unsigned len, double p, unsigned seed=DEF_SEED) {
  std::mt19937 generator;
  generator.seed(seed);

  //create the return value
  TimingStats ret;
  //initialize the random samplers
  RepeatBernoulli bern(len, p);
  BinomialShuffleOld bin_old(len, p);
  BinomialShufflePrecompute bin_new(len, p);
  PoissonOr poisson(len, p);
  FiniteDigit digit(len, p);
  BinomialShufflePrecompute bin_correction(len, digit.get_correction());
  PoissonOr poi_correction(len, digit.get_correction());
  //initialize the data storage
  std::vector<unsigned int> ber_test(n_trials);
  std::vector<unsigned int> poisson_test(n_trials);
  std::vector<unsigned int> bin_test_old(n_trials);
  std::vector<unsigned int> bin_test_new(n_trials);
  std::vector<unsigned int> digit_test(n_trials);

  //perform the test using the traditional method, store each mutation in the bernTest array
  auto begin = std::chrono::high_resolution_clock::now();
  for (unsigned i = 0; i < n_trials; ++i) {
    ber_test[i] = bern(generator);
  }
  auto end = std::chrono::high_resolution_clock::now();
  ret.bernoulli_total = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();
  
  //initialize the clock for the poisson method
  begin = std::chrono::high_resolution_clock::now();
  for (unsigned i = 0; i < n_trials; ++i) {
    poisson_test[i] = poisson(generator);
  }
  end = std::chrono::high_resolution_clock::now();
  ret.poisson_total = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();

  //initialize the clock for the old binomial shuffle
  begin = std::chrono::high_resolution_clock::now();
  for (unsigned i = 0; i < n_trials; ++i) {
    bin_test_old[i] = bin_old(generator);
  }
  end = std::chrono::high_resolution_clock::now();
  ret.binomial_old_total = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();
  
  //initialize the clock for the new binomial shuffle
  begin = std::chrono::high_resolution_clock::now();
  for (unsigned i = 0; i < n_trials; ++i) {
    bin_test_new[i] = bin_new(generator);
  }
  end = std::chrono::high_resolution_clock::now();
  ret.binomial_new_total = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();

  if (digit.get_correction() == 0) {
    begin = std::chrono::high_resolution_clock::now();
    for (unsigned i = 0; i < n_trials; ++i) {
      digit_test[i] = digit(generator);
    }
    end = std::chrono::high_resolution_clock::now();
    ret.hybrid_poisson_total = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();
    ret.hybrid_binomial_total = ret.hybrid_poisson_total;
  } else {
    //initialize the clock for the finite digit method (poisson correction)
    begin = std::chrono::high_resolution_clock::now();
    for (unsigned i = 0; i < n_trials; ++i) {
      digit_test[i] = digit(generator);
      digit_test[i] |= poi_correction(generator);
    }
    end = std::chrono::high_resolution_clock::now();
    ret.hybrid_poisson_total = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();

    //initialize the clock for the finite digit method (new binomial correction)
    begin = std::chrono::high_resolution_clock::now();
    for (unsigned i = 0; i < n_trials; ++i) {
      digit_test[i] = digit(generator);
      digit_test[i] |= bin_correction(generator);
    }
    end = std::chrono::high_resolution_clock::now();
    ret.hybrid_binomial_total = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();
  }

  std::ofstream dout;
  dout.open("dat.csv");
  //number to track how many ones there are to go in ascending order
  std::map<unsigned, OccurrenceCounter> histogram;
  for (unsigned i = 0; i < n_trials; ++i) {
    if (histogram.count(ber_test[i]) == 0) {
      histogram[ber_test[i]] = OccurrenceCounter();
    }
    if (histogram.count(poisson_test[i]) == 0) {
      histogram[poisson_test[i]] = OccurrenceCounter();
    }
    if (histogram.count(bin_test_old[i]) == 0) {
      histogram[bin_test_old[i]] = OccurrenceCounter();
    }
    if (histogram.count(bin_test_new[i]) == 0) {
      histogram[bin_test_new[i]] = OccurrenceCounter();
    }
    if (histogram.count(digit_test[i]) == 0) {
      histogram[digit_test[i]] = OccurrenceCounter();
    }
    histogram[ber_test[i]].bernoulli_occurrences += 1;
    histogram[poisson_test[i]].poisson_occurrences += 1;
    histogram[bin_test_old[i]].binomial_old_occurrences += 1;
    histogram[bin_test_new[i]].binomial_new_occurrences += 1;
    histogram[digit_test[i]].digit_occurrences += 1;
  }
  dout << "result,bernoulli,poisson,binomial old,binomial new,hybrid poisson\n";
  for (auto it = histogram.begin(); it != histogram.end(); ++it) {
    dout << it->first << ","
         << it->second.bernoulli_occurrences << ","
         << it->second.poisson_occurrences << ","
         << it->second.binomial_old_occurrences << ","
         << it->second.binomial_new_occurrences << ","
         << it->second.digit_occurrences << "," << std::endl;
  }
  dout.close();

  //calculate averages
  ret.bernoulli_avg = (double)ret.bernoulli_total / n_trials;
  ret.poisson_avg = (double)ret.poisson_total / n_trials;
  ret.binomial_old_avg = (double)ret.binomial_old_total / n_trials;
  ret.binomial_new_avg = (double)ret.binomial_new_total / n_trials;
  ret.hybrid_poisson_avg = (double)ret.hybrid_poisson_total / n_trials;
  ret.hybrid_binomial_avg = (double)ret.hybrid_binomial_total / n_trials;

  return ret;
}

int main(int argc, char** argv) {
  unsigned n_trials = NUM_TRIALS;
  unsigned len = DEF_LEN;
  unsigned seed = DEF_SEED;
  double p = PROBABILITY;
  bool silent = false;
  bool dist_test = false;

  //parse input arguments
  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i], "-p") == 0 && i != argc - 1) {
      p = atof(argv[i+1]);
      if (p < 0.0 || p > 1.0) {
        std::cout << "invalid probability" << std::endl;
        exit(0);
      }
      i++;
    }
    if (strcmp(argv[i], "-n") == 0 && i != argc - 1) {
      n_trials = atoi(argv[i+1]);
      i++;
    }
    if (strcmp(argv[i], "-l") == 0 && i != argc - 1) {
      len = atoi(argv[i+1]);
      i++;
    }
    if (strcmp(argv[i], "-s") == 0 && i != argc - 1) {
      seed = atoi(argv[i+1]);
      i++;
    }
    if (strcmp(argv[i], "-q") == 0) {
      silent = true;
    }
    if (strcmp(argv[i], "-d") == 0) {
      dist_test = true;
    }
  }

  if (!silent) {
    std::cout << "Running with:" << std::endl
              << "\ttrials:            " << n_trials << std::endl
              << "\tseed:              " << seed << std::endl
              << "\tlength:            " << len << std::endl
              << "\tprob:              " << p << std::endl;
    if (dist_test) {
      std::cout << "\ttest distribution? yes\n";
    } else {
      std::cout << "\ttest distribution? no\n";   
    }
  }

  /*unsigned int berTest[UCHAR_MAX] = {};
  unsigned int binTest[UCHAR_MAX] = {};
  unsigned int binTest2[UCHAR_MAX] = {};*/

  if (dist_test) {
    TimingStats timings = test_non_hybrid(n_trials, len, p);
    if (silent) {
      std::cout << timings.bernoulli_total  << " " << timings.bernoulli_avg
                << " " << timings.poisson_total << " " << timings.poisson_avg
                << " " << timings.binomial_old_total << " " << timings.binomial_old_avg
                << " " << timings.binomial_new_total << " " << timings.binomial_new_avg
                << " " << timings.hybrid_poisson_total << " " << timings.hybrid_poisson_avg 
                << " " << timings.hybrid_binomial_total << " " << timings.hybrid_binomial_avg << std::endl;
    } else {
      std::cout << "Bernoulli time:       " << timings.bernoulli_total
                << " \tavg: " << timings.bernoulli_avg << std::endl
                << "poisson time:         " << timings.poisson_total
                << " \tavg: " << timings.poisson_avg << std::endl
                << "binomial time (old):  " << timings.binomial_old_total
                << " \tavg: " << timings.binomial_old_avg << std::endl
                << "binomial time (new):  " << timings.binomial_new_total
                << " \tavg: " << timings.binomial_new_avg << std::endl
                << "hybrid poisson time:  " << timings.hybrid_poisson_total
                << " \tavg: " << timings.hybrid_poisson_avg << std::endl
                << "hybrid binomial time: " << timings.hybrid_binomial_total
                << " \tavg: " << timings.hybrid_binomial_avg << std::endl;
    }
  }
}
