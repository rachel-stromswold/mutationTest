#include "main.h"

//combinatorics
unsigned int choose( unsigned int n, unsigned int k ) {
  if (k > n) return 0;
  if (k * 2 > n) k = n-k;
  if (k == 0) return 1;

  unsigned int result = n;
  for( unsigned int i = 2; i <= k; ++i ) {
    result *= (n-i+1);
    result /= i;
  }
  return result;
}

//This is the function f described in the documentation.
unsigned getBitStream (unsigned n, unsigned k, unsigned x) {
  if (k == 0) { return 0; }
//  if (k == n) { return (0x01 << n)-1; }
	
  if (x < choose(n-1, k-1) ) {
    return 0x01 | (getBitStream(n-1, k-1, x) << 1);
  } else {
    return getBitStream( n-1, k, x-choose(n-1, k-1) ) << 1;
  }
}

//This is the same function f, but implemented without recursion. It is slightly faster.
unsigned getBitStream2 (unsigned n, unsigned k, unsigned x) {
  unsigned ret = 0;
  for (unsigned i = 0; k > 0 && i < n; ++i) {
/*    if (k == n-i) {
      ret = ret | ( ((0x01 << (n-i)) - 1) << i);
      break;
    }*/
    if (x < choose((n-i)-1, k-1) ) {
      ret = ret | (0x01 << i);
      k--;
    } else {
      x -= choose((n-i)-1, k-1);
    }
  }
  return ret;
}

int main(int argc, char** argv) {
  unsigned n_trials = NUM_TRIALS;
  unsigned len = DEF_LEN;
  unsigned seed = DEF_SEED;
  double p = PROBABILITY;
  bool silent = false;

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
#ifndef DIST_TEST
    if (strcmp(argv[i], "-l") == 0 && i != argc - 1) {
      len = atoi(argv[i+1]);
      i++;
    }
#endif
    if (strcmp(argv[i], "-s") == 0 && i != argc - 1) {
      seed = atoi(argv[i+1]);
      i++;
    }
    if (strcmp(argv[i], "-q") == 0) {
      silent = true;
    }
  }

  if (!silent) {
    std::cout << "Running with:" << std::endl
	      << "\ttrials: " << n_trials << std::endl
	      << "\tseed:   " << seed << std::endl
	      << "\tlength: " << len << std::endl
	      << "\tprob:   " << p << std::endl;
  }

#ifdef DIST_TEST
  unsigned int berTest[UCHAR_MAX] = {};
  unsigned int binTest[UCHAR_MAX] = {};
  unsigned int binTest2[UCHAR_MAX] = {};
#endif
  std::mt19937 generator;
  generator.seed(seed);

  //perform the test using the traditional method, store each mutation in the bernTest array
  auto begin = std::chrono::high_resolution_clock::now();
  std::bernoulli_distribution bern(p);
  for (unsigned i = 0; i < n_trials; ++i) {
    unsigned berVal = 0;
    for (unsigned j = 0; j < len; ++j) {
      if ( bern( generator ) ) {
	berVal = berVal | (0x01 << j);
      }
    }
#ifdef DIST_TEST
    berTest[berVal]++;
#endif
  }
  auto end = std::chrono::high_resolution_clock::now();
  auto berDuration = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();

  //initialize the clock for the new method with recursion
  begin = std::chrono::high_resolution_clock::now();
  //use the modified procedure, store each mutation mask in the binTest array
  std::binomial_distribution<unsigned> bin(len, p);
  for (unsigned i = 0; i < n_trials; ++i) {
    unsigned num_ones = bin(generator);

    //minus 1 because we start indexing from 0
    std::uniform_int_distribution<unsigned> dist(0, choose(len, num_ones) - 1);

#ifdef DIST_TEST
    unsigned binVal = getBitStream(len, num_ones, dist(generator));
    binTest[binVal]++;
#else
    getBitStream(len, num_ones, dist(generator));
#endif
  }
  end = std::chrono::high_resolution_clock::now();
  auto binDuration = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();
  
  //initialize the clock for the new method without recursion
  begin = std::chrono::high_resolution_clock::now();
  for (unsigned i = 0; i < n_trials; ++i) {
    unsigned num_ones = bin(generator);

    std::uniform_int_distribution<unsigned char> dist(0, choose(len, num_ones) - 1);
#ifdef DIST_TEST
    unsigned binVal = getBitStream2(len, num_ones, dist(generator));
    binTest2[binVal]++;
#else
    getBitStream2(len, num_ones, dist(generator));
#endif
  }
  end = std::chrono::high_resolution_clock::now();
  auto binDuration2 = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();

#ifdef DIST_TEST
  std::ofstream dout;
  dout.open("dat.csv");
  //number to track how many ones there are to go in ascending order
  unsigned m = 0;
  unsigned lastShift = 0;
  unsigned nextShift = 1;
  for (unsigned i = 0; i <= UCHAR_MAX; ++i) {
    //only print data out if there actually were results
    if (binTest[i] || berTest[i]) {
      if (i < 10) {
	std::cout << " ";
      }
      if (i < 100) {
	std::cout << " ";
      }
      std::cout << i << " Bernoulli";
      unsigned j = 0;
      //print results for the bernoulli tests
      for (; j < berTest[i]*N_STARS/n_trials; ++j) {
	std::cout << "*";
      }
      //print results for the binomial tests
      std::cout << std::endl << "    Binomial ";
      for (unsigned j = 0; j < binTest[i]*N_STARS/n_trials; ++j) {
	std::cout << "*";
      }
      //print results for the other binomial tests
      std::cout << std::endl << "    Binomial2";
      for (unsigned j = 0; j < binTest[i]*N_STARS/n_trials; ++j) {
	std::cout << "*";
      }
      std::cout << std::endl;
    }
    if (i == nextShift) {
      m++;
      lastShift = nextShift;
      nextShift += choose(8, m);
    }
    //update the outputted distribution graph
    //masks are tracked in descending order of the number of 1 bits
    int ind = getBitStream(8, m, i-lastShift);
    dout << ind << "," << berTest[ind] << "," << binTest[ind] << "," << binTest2[ind] << std::endl;
  }
  dout.close();
#endif //DIST_TEST
  if (silent) {
    std::cout << " " << berDuration << " " << double(berDuration)/n_trials
	      << " " << binDuration << " " << double(binDuration)/n_trials
	      << " " << binDuration2 << " " << double(binDuration2)/n_trials << " " << std::endl;
  } else {
    std::cout << "Bernoulli time              : " << berDuration
	      << "\tavg: " << double(berDuration)/n_trials << std::endl
	      << "binomial time (w/ recursion): " << binDuration
	      << "\tavg: " << double(binDuration)/n_trials << std::endl
	      << "binomial time (wo recursion): " << binDuration2
	      << "\tavg: " << double(binDuration2)/n_trials << std::endl;
  }
}
