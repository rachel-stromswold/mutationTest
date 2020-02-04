#include "main.h"

#define TEST_DIST 0
#define TEST_PERC 1
#define TEST_PROB 2

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
  unsigned hybrid_poisson_total = 0;
  unsigned hybrid_binomial_total = 0;

  double bernoulli_avg = 0;
  double poisson_avg = 0;
  double binomial_old_avg = 0;
  double binomial_new_avg = 0;
  double hybrid_poisson_avg = 0;
  double hybrid_binomial_avg = 0;

  void calc_averages(_uint n_trials) {
    bernoulli_avg = (double)bernoulli_total / n_trials;
    poisson_avg = (double)poisson_total / n_trials;
    binomial_old_avg = (double)binomial_old_total / n_trials;
    binomial_new_avg = (double)binomial_new_total / n_trials;
    hybrid_poisson_avg = (double)hybrid_poisson_total / n_trials;
    hybrid_binomial_avg = (double)hybrid_binomial_total / n_trials;
  }
};

TimingStats test_non_hybrid(unsigned n_trials, unsigned len, double p, unsigned seed=DEF_SEED) {
  std::mt19937 generator;
  generator.seed(seed);

  //create the return value
  TimingStats ret;
  //initialize the random samplers and the data storage
  SampleData<RepeatBernoulli> test_bern(len, p, n_trials);
    TimingTracker<SampleData<RepeatBernoulli>, std::mt19937> track_bern(n_trials);
  SampleData<PoissonOr> test_pois(len, p, n_trials);
    TimingTracker<SampleData<PoissonOr>, std::mt19937> track_pois(n_trials);
  SampleData<BinomialShuffleOld> test_bin_old(len, p, n_trials);
    TimingTracker<SampleData<BinomialShuffleOld>, std::mt19937> track_bin_old(n_trials);
  SampleData<BinomialShufflePrecompute> test_bin_new(len, p, n_trials);
    TimingTracker<SampleData<BinomialShufflePrecompute>, std::mt19937> track_bin_new(n_trials);
  SampleData<FiniteDigit<PoissonOr>> test_hyp(len, p, n_trials);
    TimingTracker<SampleData<FiniteDigit<PoissonOr>>, std::mt19937> track_hyp(n_trials);
  SampleData<FiniteDigit<BinomialShufflePrecompute>> test_hyb(len, p, n_trials);
    TimingTracker<SampleData<FiniteDigit<BinomialShufflePrecompute>>, std::mt19937> track_hyb(n_trials);

  std::cout << "Now running bernoulli test\n";
  ret.bernoulli_total = track_bern.run(test_bern, generator);
  std::cout << "Now running poisson test\n";
  ret.poisson_total = track_pois.run(test_pois, generator);
  std::cout << "Now running old-binomial test\n";
  ret.binomial_old_total = track_bin_old.run(test_bin_old, generator);
  std::cout << "Now running new-binomial test\n";
  ret.binomial_new_total = track_bin_new.run(test_bin_new, generator); 
  std::cout << "Now running poisson finite digit test\n";
  ret.hybrid_poisson_total = track_hyp.run(test_hyp, generator); 
  std::cout << "Now running binomial finite digit test\n";
  ret.hybrid_binomial_total = track_hyb.run(test_hyb, generator); 

  std::ofstream dout;
  dout.open("dat.csv");
  //number to track how many ones there are to go in ascending order
  std::map<unsigned, OccurrenceCounter> histogram;
  for (unsigned i = 0; i < n_trials; ++i) {
    if (histogram.count(test_bern.get(i)) == 0) {
      histogram[test_bern.get(i)] = OccurrenceCounter();
    }
    if (histogram.count(test_pois.get(i)) == 0) {
      histogram[test_pois.get(i)] = OccurrenceCounter();
    }
    if (histogram.count(test_bin_old.get(i)) == 0) {
      histogram[test_bin_old.get(i)] = OccurrenceCounter();
    }
    if (histogram.count(test_bin_new.get(i)) == 0) {
      histogram[test_bin_new.get(i)] = OccurrenceCounter();
    }
    if (histogram.count(test_hyb.get(i)) == 0) {
      histogram[test_hyb.get(i)] = OccurrenceCounter();
    }
    histogram[test_bern.get(i)].bernoulli_occurrences += 1;
    histogram[test_pois.get(i)].poisson_occurrences += 1;
    histogram[test_bin_old.get(i)].binomial_old_occurrences += 1;
    histogram[test_bin_new.get(i)].binomial_new_occurrences += 1;
    histogram[test_hyb.get(i)].digit_occurrences += 1;
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
  ret.calc_averages(n_trials);
  return ret;
}

typedef PercolationTracker<RepeatBernoulli> PercBern;
typedef PercolationTracker<PoissonOr> PercPois;
typedef PercolationTracker<BinomialShuffleOld> PercBinOld;
typedef PercolationTracker<BinomialShufflePrecompute> PercBinNew;
typedef PercolationTracker<FiniteDigit<PoissonOr>> PercHyp;
typedef PercolationTracker<FiniteDigit<BinomialShufflePrecompute>> PercHyb;

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

  PercBern perc_bern(n_trials, p, relax, len, t_max);
    TimingTracker<PercBern, std::mt19937> track_bern(t_max);
  PercPois perc_pois(n_trials, p, relax, len, t_max);
    TimingTracker<PercPois, std::mt19937> track_pois(t_max);
  PercBinOld perc_bin_old(n_trials, p, relax, len, t_max);
    TimingTracker<PercBinOld, std::mt19937> track_bin_old(t_max);
  PercBinNew perc_bin_new(n_trials, p, relax, len, t_max);
    TimingTracker<PercBinNew, std::mt19937> track_bin_new(t_max);
  PercHyp perc_hyp(n_trials, p, relax, len, t_max);
    TimingTracker<PercHyp, std::mt19937> track_hyp(t_max);
  PercHyb perc_hyb(n_trials, p, relax, len, t_max);
    TimingTracker<PercHyb, std::mt19937> track_hyb(t_max);

  
  std::cout << "Now running bernoulli percolation\n";
  ret.bernoulli_total = track_bern.run(perc_bern, generator);
  std::cout << "Now running poisson percolation\n";
  ret.poisson_total = track_pois.run(perc_pois, generator);
  std::cout << "Now running old-binomial percolation\n";
  ret.binomial_old_total = track_bin_old.run(perc_bin_old, generator);
  std::cout << "Now running new-binomial percolation\n";
  ret.binomial_new_total = track_bin_new.run(perc_bin_new, generator); 
  std::cout << "Now running poisson finite digit percolation\n";
  ret.hybrid_poisson_total = track_hyp.run(perc_hyp, generator); 
  std::cout << "Now running binomial finite digit percolation\n";
  ret.hybrid_binomial_total = track_hyb.run(perc_hyb, generator);

  ret.calc_averages(n_trials);

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

TimingStats test_probabilities(unsigned n_trials, unsigned len, double p, unsigned seed=DEF_SEED) {
  std::mt19937 generator;
  generator.seed(seed);
  TimingStats ret;

  BitsData<RepeatBernoulli> bern(len, p);
    TimingTracker<BitsData<RepeatBernoulli>, std::mt19937> track_bern(n_trials);
  BitsData<PoissonOr> pois(len, p);
    TimingTracker<BitsData<PoissonOr>, std::mt19937> track_pois(n_trials);
  BitsData<BinomialShuffleOld> bin_old(len, p);
    TimingTracker<BitsData<BinomialShuffleOld>, std::mt19937> track_bin_old(n_trials);
  BitsData<BinomialShufflePrecompute> bin_new(len, p);
    TimingTracker<BitsData<BinomialShufflePrecompute>, std::mt19937> track_bin_new(n_trials);
  BitsData<FiniteDigit<PoissonOr>> hyp(len, p);
    TimingTracker<BitsData<FiniteDigit<PoissonOr>>, std::mt19937> track_hyp(n_trials);
  BitsData<FiniteDigit<BinomialShufflePrecompute>> hyb(len, p);
    TimingTracker<BitsData<FiniteDigit<BinomialShufflePrecompute>>, std::mt19937> track_hyb(n_trials);

  std::cout << "Now running bernoulli probabilities\n";
  ret.bernoulli_total = track_bern.run(bern, generator);
  std::cout << "Now running poisson probabilities\n";
  ret.poisson_total = track_pois.run(pois, generator);
  std::cout << "Now running old-binomial probabilities\n";
  ret.binomial_old_total = track_bin_old.run(bin_old, generator);
  std::cout << "Now running new-binomial probabilities\n";
  ret.binomial_new_total = track_bin_new.run(bin_new, generator); 
  std::cout << "Now running poisson finite digit probabilities\n";
  ret.hybrid_poisson_total = track_hyp.run(hyp, generator); 
  std::cout << "Now running binomial finite digit probabilities\n";
  ret.hybrid_binomial_total = track_hyb.run(hyb, generator);

  std::cout << "bernoulli:\n";
  bern.print_statistics(p);
  std::cout << "poisson:\n";
  pois.print_statistics(p);
  std::cout << "binomial old:\n";
  bin_old.print_statistics(p);
  std::cout << "binomial new:\n";
  bin_new.print_statistics(p);
  std::cout << "finite digit (new binomial):\n";
  hyb.print_statistics(p);
  std::cout << "finite digit (poisson):\n";
  hyp.print_statistics(p);

  ret.calc_averages(n_trials);
  return ret;
}

int main(int argc, char** argv) {
  unsigned n_trials = NUM_TRIALS;
  unsigned n_steps = NUM_STEPS;
  unsigned len = DEF_LEN;
  unsigned seed = DEF_SEED;
  double p = PROBABILITY;
  bool silent = false;
  _uint test_type = TEST_PERC;
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
    if (strcmp(argv[i], "--distribution") == 0) {
      test_type = TEST_DIST;
    }
    if (strcmp(argv[i], "--probabilities") == 0) {
      test_type = TEST_PROB;
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
    if (test_type == TEST_DIST) {
      std::cout << "\ttest type:         distribution\n";
    } else if (test_type == TEST_PROB) {
      std::cout << "\ttest type:         probability\n";
    } else {
      std::cout << "\ttest type:         percolation\n";
    }
  }

  /*unsigned int berTest[UCHAR_MAX] = {};
  unsigned int binTest[UCHAR_MAX] = {};
  unsigned int binTest2[UCHAR_MAX] = {};*/

  TimingStats timings;
  if (test_type == TEST_DIST) {
    timings = test_non_hybrid(n_trials, len, p);
  } else if (test_type == TEST_PERC) {
    timings = run_percolation(n_trials, len, p, relax, n_steps, seed);
  } else if (test_type == TEST_PROB) {
    timings = test_probabilities(n_trials, len, p, seed);
  }
  if (silent) {
    std::cout << timings.bernoulli_total  << " " << timings.bernoulli_avg
	      << " " << timings.poisson_total << " " << timings.poisson_avg
	      << " " << timings.binomial_old_total << " " << timings.binomial_old_avg
	      << " " << timings.binomial_new_total << " " << timings.binomial_new_avg
	      << " " << timings.hybrid_poisson_total << " " << timings.hybrid_poisson_avg 
	      << " " << timings.hybrid_binomial_total << " " << timings.hybrid_binomial_avg << std::endl;
  } else {
    std::cout << "Bernoulli time:       " << int_fmt(timings.bernoulli_total, D_LEN)
	      << " \tavg: " << int_fmt(timings.bernoulli_avg, D_LEN) << std::endl
	      << "poisson time:         " << int_fmt(timings.poisson_total, D_LEN)
	      << " \tavg: " << int_fmt(timings.poisson_avg, D_LEN) << std::endl
	      << "binomial time (old):  " << int_fmt(timings.binomial_old_total, D_LEN)
	      << " \tavg: " << int_fmt(timings.binomial_old_avg, D_LEN) << std::endl
	      << "binomial time (new):  " << int_fmt(timings.binomial_new_total, D_LEN)
	      << " \tavg: " << int_fmt(timings.binomial_new_avg, D_LEN) << std::endl
	      << "hybrid poisson time:  " << int_fmt(timings.hybrid_poisson_total, D_LEN)
	      << " \tavg: " << int_fmt(timings.hybrid_poisson_avg, D_LEN) << std::endl
	      << "hybrid binomial time: " << int_fmt(timings.hybrid_binomial_total, D_LEN)
	      << " \tavg: " << int_fmt(timings.hybrid_binomial_avg, D_LEN) << std::endl;
  }
}
