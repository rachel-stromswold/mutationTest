#include "main.h"

//namespace plt = matplotlibcpp;

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
  FiniteDigit<BinomialShufflePrecompute> digit_bin(len, p);
  FiniteDigit<PoissonOr> digit_poi(len, p);
  /*BinomialShufflePrecompute bin_correction(len, digit.get_correction());
  PoissonOr poi_correction(len, digit.get_correction());*/
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

  //initialize the clock for the finite digit method (poisson correction)
  begin = std::chrono::high_resolution_clock::now();
  for (unsigned i = 0; i < n_trials; ++i) {
    digit_test[i] = digit_poi(generator);
    //digit_test[i] |= poi_correction(generator);
  }
  end = std::chrono::high_resolution_clock::now();
  ret.hybrid_poisson_total = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();

  //initialize the clock for the finite digit method (new binomial correction)
  begin = std::chrono::high_resolution_clock::now();
  for (unsigned i = 0; i < n_trials; ++i) {
    digit_test[i] = digit_bin(generator);
    //digit_test[i] |= bin_correction(generator);
  }
  end = std::chrono::high_resolution_clock::now();
  ret.hybrid_binomial_total = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();
 

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

TimingStats run_percolation(unsigned n_trials, unsigned len, double p, bool relax, unsigned t_max=10000, unsigned seed=DEF_SEED) {
  std::mt19937 generator;
  generator.seed(seed);

  TimingStats ret;

  unsigned bernoulli_total = 0;
  unsigned poisson_total = 0;
  unsigned binomial_old_total = 0;
  unsigned binomial_new_total = 0;
  double hybrid_poisson_total = 0;
  double hybrid_binomial_total = 0;

  PercolationTracker<RepeatBernoulli> perc_bern(n_trials, p, relax, len, t_max);
  PercolationTracker<PoissonOr> perc_pois(n_trials, p, relax, len, t_max);
  PercolationTracker<BinomialShuffleOld> perc_bin_old(n_trials, p, relax, len, t_max);
  PercolationTracker<BinomialShufflePrecompute> perc_bin_new(n_trials, p, relax, len, t_max);
  PercolationTracker<FiniteDigit<PoissonOr>> perc_hyp(n_trials, p, relax, len, t_max);
  PercolationTracker<FiniteDigit<BinomialShufflePrecompute>> perc_hyb(n_trials, p, relax, len, t_max);

  _uint interval = 500;
  std::cout << "Now running bernoulli test\nv";
  for (_uint i = 0; i < t_max/INTERVAL; ++i) {
    std::cout << " ";
  }
  std::cout << "v\n ";
  auto begin = std::chrono::high_resolution_clock::now();
  for (_uint i = 0; i < t_max; ++i) {
    perc_bern.update(generator);
    if (i % INTERVAL == 0) {
      std::cout << "#" << std::flush;
    }
  }
  auto end = std::chrono::high_resolution_clock::now();
  std::cout << "perc bern\n\t";perc_bern.print_probabilities();
  ret.bernoulli_total = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();
  std::cout << "\n";

  std::cout << "Now running poisson test\nv";
  for (_uint i = 0; i < t_max/INTERVAL; ++i) {
    std::cout << " ";
  }
  std::cout << "v\n ";
  begin = std::chrono::high_resolution_clock::now();
  for (_uint i = 0; i < t_max; ++i) {
    perc_pois.update(generator);
    if (i % INTERVAL == 0) {
      std::cout << "#" << std::flush;
    }
  }
  end = std::chrono::high_resolution_clock::now();
  std::cout << "perc pois\n\t";perc_pois.print_probabilities();
  ret.poisson_total = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();
  std::cout << "\n";

  std::cout << "Now running old-binomial test\nv";
  for (_uint i = 0; i < t_max/INTERVAL; ++i) {
    std::cout << " ";
  }
  std::cout << "v\n ";
  begin = std::chrono::high_resolution_clock::now();
  for (_uint i = 0; i < t_max; ++i) {
    perc_bin_old.update(generator);
    if (i % INTERVAL == 0) {
      std::cout << "#" << std::flush;
    }
  }
  end = std::chrono::high_resolution_clock::now();
  std::cout << "perc bino\n\t";perc_bin_old.print_probabilities();
  ret.binomial_old_total = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();
  std::cout << "\n";

  std::cout << "Now running new-binomial test\nv";
  for (_uint i = 0; i < t_max/INTERVAL; ++i) {
    std::cout << " ";
  }
  std::cout << "v\n ";
  begin = std::chrono::high_resolution_clock::now();
  for (_uint i = 0; i < t_max; ++i) {
    perc_bin_new.update(generator);
    if (i % INTERVAL == 0) {
      std::cout << "#" << std::flush;
    }
  }
  end = std::chrono::high_resolution_clock::now();
  std::cout << "perc binn\n\t";perc_bin_new.print_probabilities();
  ret.binomial_new_total = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();
  std::cout << "\n";

  std::cout << "Now running poisson finite digit test\nv";
  for (_uint i = 0; i < t_max/INTERVAL; ++i) {
    std::cout << " ";
  }
  std::cout << "v\n ";
  begin = std::chrono::high_resolution_clock::now();
  for (_uint i = 0; i < t_max; ++i) {
    perc_hyp.update(generator);
    if (i % INTERVAL == 0) {
      std::cout << "#" << std::flush;
    }
  }
  end = std::chrono::high_resolution_clock::now();
  std::cout << "perc hyp\n\t";perc_hyp.print_probabilities();
  ret.hybrid_poisson_total = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();
  std::cout << "\n";

  std::cout << "Now running binomial finite digit test\nv";
  for (_uint i = 0; i < t_max/INTERVAL; ++i) {
    std::cout << " ";
  }
  std::cout << "v\n ";
  begin = std::chrono::high_resolution_clock::now();
  for (_uint i = 0; i < t_max; ++i) {
    perc_hyb.update(generator);
    if (i % INTERVAL == 0) {
      std::cout << "#" << std::flush;
    }
  }
  end = std::chrono::high_resolution_clock::now();
  std::cout << "perc hyb\n\t";perc_hyb.print_probabilities();
  ret.hybrid_binomial_total = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();
  std::cout << "\n";

  ret.bernoulli_avg = (double)ret.bernoulli_total / n_trials;
  ret.poisson_avg = (double)ret.poisson_total / n_trials;
  ret.binomial_old_avg = (double)ret.binomial_old_total / n_trials;
  ret.binomial_new_avg = (double)ret.binomial_new_total / n_trials;
  ret.hybrid_poisson_avg = (double)ret.hybrid_poisson_total / n_trials;
  ret.hybrid_binomial_avg = (double)ret.hybrid_binomial_total / n_trials;

  std::vector<double> time = perc_bern.get_time_arr();
  _uint n = time.size();

  double* ln_time = new double[n];
  double* ln_occ = new double[n];
  double* ln_rho = new double[n];
  std::vector<double> occ_fitted(n);
  std::vector<double> rho_fitted(n);

  double theta_fit, theta_cov, theta_sumsq_res, theta_r2, occ_mean, occ_ss;
  double beta_fit, beta_cov, beta_sumsq_res, beta_r2, rho_mean, rho_ss;
  std::ofstream dout;
  dout.open("occupations.csv");
  dout << "time,n(t),rho(t),n_surviving\n";
  for (_uint t = 0; t < n; ++t) {
    double avg_n = ( perc_bern.get_n(t) + perc_pois.get_n(t) + perc_bin_old.get_n(t) + perc_bin_new.get_n(t) + perc_hyp.get_n(t) + perc_hyb.get_n(t) ) / 6;
    double avg_rho = ( perc_bern.get_rho(t) + perc_pois.get_rho(t) + perc_bin_old.get_rho(t) + perc_bin_new.get_rho(t) + perc_hyp.get_rho(t) + perc_hyb.get_rho(t) ) / 6;
    double avg_survivors = ( perc_bern.get_survivors(t) + perc_pois.get_survivors(t) + perc_bin_old.get_survivors(t) + perc_bin_new.get_survivors(t) + perc_hyp.get_survivors(t) + perc_hyb.get_survivors(t) ) / 6;

    ln_time[t] = (time[t] <= 0)? 0 : log(time[t]);
    ln_occ[t] = log(avg_n);
    ln_rho[t] = log(avg_rho);
    occ_mean += ln_occ[t] / n;
    rho_mean += ln_rho[t] / n;
    dout << t << ","
         << avg_n << "," << avg_rho << "," << avg_survivors << ","
         << perc_bern.get_n(t) << "," << perc_bern.get_rho(t) << "," << perc_bern.get_survivors(t) << ","
	 << perc_pois.get_n(t) << "," << perc_pois.get_rho(t) << "," << perc_pois.get_survivors(t) << ","
	 << perc_bin_old.get_n(t) << "," << perc_bin_old.get_rho(t) << "," << perc_bin_old.get_survivors(t) << ","
	 << perc_bin_new.get_n(t) << "," << perc_bin_new.get_rho(t) << "," << perc_bin_new.get_survivors(t) << ","
	 << perc_hyp.get_n(t) << "," << perc_hyp.get_rho(t) << "," << perc_hyp.get_survivors(t) << ","
	 << perc_hyb.get_n(t) << "," << perc_hyb.get_rho(t) << "," << perc_hyb.get_survivors(t) << ","
	 << std::endl;
  }

  gsl_fit_mul(ln_time, 1, ln_occ, 1, time.size(), &theta_fit, &theta_cov, &theta_sumsq_res);
  gsl_fit_mul(ln_time, 1, ln_rho, 1, time.size(), &beta_fit, &beta_cov, &beta_sumsq_res);

  for (_uint i = 0; i < n; ++i) {
    occ_ss += pow(ln_occ[i] - occ_mean, 2);
    rho_ss += pow(ln_rho[i] - rho_mean, 2);
    occ_fitted[i] = pow(time[i], theta_fit);
    rho_fitted[i] = pow(time[i], beta_fit);
  }
  theta_r2 = 1.0 - theta_sumsq_res/occ_ss;
  beta_r2 = 1.0 - beta_sumsq_res/rho_ss;
  double theta_r2_adj = (1.0-(1.0-theta_r2)*(n-1)/(n-2));
  double beta_r2_adj = (1.0-(1.0-beta_r2)*(n-1)/(n-2));

  std::cout << "data\t| power law\t| std_dev\t| R^2\t|R^2 (adjusted)\n"
	    << "n(t)\t| " << theta_fit << "\t| " << sqrt(theta_cov) << "\t| " << theta_r2 << "\t| " << theta_r2_adj << "\n"
	    << "rho(t)\t| " << beta_fit << "\t| " << sqrt(beta_cov) << "\t| " << beta_r2 << "\t| " << beta_r2_adj << "\n";

  delete[] ln_time;
  delete[] ln_occ;
  delete[] ln_rho;

  return ret;
}

int main(int argc, char** argv) {
  unsigned n_trials = NUM_TRIALS;
  unsigned n_steps = NUM_STEPS;
  unsigned len = DEF_LEN;
  unsigned seed = DEF_SEED;
  double p = PROBABILITY;
  bool silent = false;
  bool dist_test = false;
  bool relax = false;

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
    if (strcmp(argv[i], "-t") == 0 && i != argc - 1) {
      n_steps = atoi(argv[i+1]);
      i++;
    }
    if (strcmp(argv[i], "-q") == 0) {
      silent = true;
    }
    if (strcmp(argv[i], "-d") == 0) {
      dist_test = true;
    }
    if (strcmp(argv[i], "-r") == 0 && i != argc - 1) {
      relax = true;
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

  TimingStats timings;
  if (dist_test) {
    timings = test_non_hybrid(n_trials, len, p);
  } else {
    timings = run_percolation(n_trials, len, p, relax, n_steps, seed);
  }
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
