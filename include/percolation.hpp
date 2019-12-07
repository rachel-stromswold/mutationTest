#include "matplotlibcpp.h"
#include "bitstrings.hpp"

template <typename Sampler, unsigned STRING_SIZE>
class DirectedPercolator1() {
private:
  const _uint last_bit = 1 << (STRING_SIZE - 1);
  _uint t_max;
  _uint t = 0;
  _uint steady_state_t = t_max - 1;
  Sampler s;
  double p;
  std::vector<std::vector<_uint>> bitstrings;//bitstrings[time][spatial location]

  _uint one_state(_uint t = STRING_SIZE) {
    if (t >= STRING_SIZE) { return (last_bit-1) | last_bit; }
    return (1 << t) - 1;
  }

public:
  DirectedPercolator1(_uint pt_max, double p_p) {
    _uint n_strings = 0;
    for (_uint i = 0; i < t_max; ++i) {
      if (i % STRING_SIZE == 0) {
        ++n_strings;
      }
      bitstrings[i].resize(n_strings);
    }
    bitstrings[0] = 1;
  }

  template <typename Generator>
  void update(Generator& g) {
    ++t;
    if (t < steady_state_t && t < t_max) {
      _uint l_bonds = s(g);
      _uint r_bonds = s(g);
      bitstrings[t][i] = (bitstrings[t-1][0] & r_bonds) | ( (bitstrings[t-1][0] & l_bonds) << 1 );

      bool all_zeros = (bitstrings[t][0] == 0);
      bool all_ones = (bitstrings[t][0] == one_state(t));

      for (_uint i = 1; i < bitstrings[t].size(); ++i) {
        //for the time being we need to keep track of the left bonds to handle carryover
        _uint new_l_bonds = s(g);
        r_bonds = s(g);
        //handle carryover
        bitstrings[t][i] = (bitstrings[t-1][i-1] & (l_bonds & last_bit)) >> (STRING_SIZE - 1);
        l_bonds = new_l_bonds;
        if (i < bitstrings[t-1].size) {
          bitstrings[t][i] |= (bitstrings[t-1][i] & r_bonds) | ( (bitstrings[t-1][i] & l_bonds) << 1 );
        }
        if (bitstrings[t][i] != 0) { all_zeros = false; }
        if (i < bitstrings[t].size()-1 && bitstrings[t][i] != one_state(t)) { all_ones = false; }
      }
      if (bitstrings[t][bitstrings[t].size() - 1] < one_state(t % STRING_SIZE)) { all_ones = false; }
      if (all_zeros || all_ones) {
        steady_state_t = t;
      }
    }
  }

  std::vector<double> get_n_t() {
    std::vector<_uint> result(t_max);
    for (_uint tt = 0; tt <= steady_state_t; ++tt) {
      _uint count = 0;
      for (_uint i = 0; i < bitstrings[t].size(); ++i) {
        _uint tmp = bitstrings[t][i]
        while (tmp > 0) {
          count += tmp & 1;
          tmp = tmp >> 1;
        }
      }
      result[tt] = count;
    }
    for (_uint tt = steady_state_t + 1; tt < t_max; ++tt) {
      result[tt] = result[steady_state_t];
    }
    return result;
  }

  std::vector<double> get_rho_t() {
    std::vector<double> counts = get_n_t();
    std::vector<double> result(t_max);
    for (_uint tt = 0; tt < result.size(); ++tt) {
      result[tt] = (double)counts[tt] / (tt + 1);
    }
    return result;
  }
};

template <typename Sampler, unsigned STRING_SIZE>
class PercolationTracker {
private:
  std::vector<_uint> avg_occupation_n;
  std::vector<double> avg_rho;
  Sampler s;
  std::vector<DirectedPercolator1<Sampler, STRING_SIZE>> percs;

public:
  PercolationTracker(_uint n_samples);
  void update();
};
