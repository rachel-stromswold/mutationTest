#ifndef MAIN_H
#define MAIN_H

//#include "include/bitstrings.hpp"
#include "include/percolation.hpp"
#include "include/probability.hpp"
//#include <gsl/gsl_fit.h>
#include <map>
#include <iostream>
#include <string.h>
#include <fstream>
//#include <climits>
#include <chrono>

#define NUM_TRIALS	1000000
#define NUM_STEPS	32768
#define PROBABILITY	0.01
#define COL_WIDTH	80
#define STATUS_BAR_LEN	40
#define D_LEN 10

#define N_STARS		200

#define DEF_LEN		32
#define DEF_SEED	3141593

std::string int_fmt(long int val, unsigned char len) {
  if (len < 5) { len = 5; }
  std::string ret(len+1, ' ');
  int n_digits = (int)ceil(log(val) / log(10));
  
  if (n_digits > len) {
    ret[len-1] = n_digits % 10;
    ret[len-2] = n_digits / 10;
    ret[len-3] = '+';
    ret[len-4] = 'e';
    ret[0] = val / ( (int)pow(10, n_digits-1) );
    if (len > 5) { ret[1] = '.'; }
    if (len > 6) {
      for (_uint i = 0; i < len - 6; ++i) {
        ret[len - i - 5] = (unsigned char)(val % 10) + '0';
        val /= 10;
      }
    }
  } else {
    for (_uint i = len; i > len - n_digits; --i) {
      ret[i] = (unsigned char)(val % 10) + '0';
      val /= 10;
    }
  }
  return ret;
} 

//takes a class Ticker which has a function called update(U) 
template <typename Ticker, typename U>
class TimingTracker {
private:
  _uint n_steps;
  _uint interval;

public:
  TimingTracker(_uint pn_steps, _uint p_interval = 0) : n_steps(pn_steps), interval(p_interval) {
    if (interval == 0) {
      interval = (n_steps + STATUS_BAR_LEN-1) / STATUS_BAR_LEN;
    }
  }

  unsigned run(Ticker& tick, U& u) {
    std::cout << "v";
    for (_uint i = 0; i < n_steps / interval; ++i) {
      std::cout << " ";
    }
    if (n_steps % interval != 0) {
      std::cout << " ";
    }
    std::cout << "v\n ";

    //initialize the clock for the new binomial shuffle
    auto begin = std::chrono::high_resolution_clock::now();
    //step through
    for (_uint i = 0; i < n_steps; ++i) {
      tick.update(u);
      if (i % interval == 0) {
        std::cout << "#" << std::flush;
      }
    }
    //finalize the clock
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "\n";
    //return the total time taken
    return std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();
  }
};

template<typename Sampler>
class SampleData {
private:
  std::vector<unsigned int> data;
  Sampler s;
  _uint i = 0;
public:
  SampleData(_uint p_len, double p_p, _uint pn_trials) : s(p_len, p_p),  data(pn_trials) {}

  void update(std::mt19937& generator) {
    if (i >= data.size()) {
      data.push_back(s(generator));
    } else {
      data[i] = s(generator);
    }
    ++i;
  }

  unsigned int get(_uint ind) { return data[ind]; }
};

#endif //MAIN_H
