#ifndef BITSTRINGS_H
#define BITSTRINGS_H

#include <random>
#include <iostream>
#include <math.h>

typedef size_t _uint;

//combinatorics
_uint choose( _uint n, _uint k ) {
  if (k > n) return 0;
  if (k * 2 > n) k = n-k;
  if (k == 0) return 1;

  _uint result = n;
  for( _uint i = 2; i <= k; ++i ) {
    result *= (n-i+1);
    result /= i;
  }
  return result;
}

//This is the function f described in the documentation.
_uint get_bit_stream_slow(_uint n, _uint k, _uint x) {
  if (k == 0) { return 0; }
//  if (k == n) { return (0x01 << n)-1; }
	
  if (x < choose(n-1, k-1) ) {
    return 0x01 | (get_bit_stream_slow(n-1, k-1, x) << 1);
  } else {
    return get_bit_stream_slow( n-1, k, x-choose(n-1, k-1) ) << 1;
  }
}

class RepeatBernoulli {
private:
  _uint n_;
  double p_;
  std::bernoulli_distribution bern;

public:
  RepeatBernoulli(_uint n, double p) : bern(p) {
    n_ = n;
    p_ = p;
  }

  void reset() { bern.reset(); }
  _uint n() { return n_; }
  double p() { return p_; }

  template <class Generator>
  _uint operator()(Generator& g) {
    _uint berVal = 0;
    for (_uint j = 0; j < n_; ++j) {
      if ( bern(g) ) {
        berVal = berVal | (0x01 << j);
      }
    }
    
    return berVal;
  }
};

class FiniteDigit {
private:
  _uint n_;
  double p_;
  std::vector<bool> binary_rep;

public:
  FiniteDigit(_uint n, double p, _uint accuracy=3) : binary_rep(accuracy) {
    n_ = n;
    p_ = p;
    for (_uint i = 0; i < binary_rep.size(); ++i) {
      if (p >= 0.5) {
        binary_rep[i] = 1;
        p = p*2 - 1;
      } else {
        binary_rep[i] = 0;
        p = p*2;
      }
    }
    while (binary_rep.size() > 0 && binary_rep[binary_rep.size()-1] == 0) {
      binary_rep.pop_back();
    }
  }
  double get_correction() {
    double actual_p = 0;
    for (_uint i = 0; i < binary_rep.size(); ++i) {
      if (binary_rep[i]) {
        actual_p += 0.5/((double)(1 << i));
      }
    }
    if (actual_p == 1) { return 0; }
    return (p_ - actual_p)/(1 - actual_p);
  }
  
  template <class Generator>
  _uint operator()(Generator& g) {
    _uint val = 0;

    std::uniform_int_distribution<_uint> unif(0, (1 << n_)-1);
    //std::cout << val;
    for (int i = binary_rep.size()-1; i >= 0; --i) {
      if (binary_rep[i]) {
        val = val | unif(g);
        //std::cout << " ored " << val;
      } else {
        val = val & unif(g);
        //std::cout << " anded " << val;
      }
    }
    //std::cout << std::endl;
    return val;
  }
};

//implementation found in Watanabe et al.
class BinomialShuffleOld {
private:
  _uint n_;
  double p_;
  std::binomial_distribution<_uint> bin;
  _uint invert = 0;

public:
  BinomialShuffleOld(_uint n, double p) : bin(n, p) {
    n_ = n;
    if (p > 0.5) {
      p_ = 1 - p;
      invert = (1 << n) - 1;
      bin = std::binomial_distribution<_uint>(n, p_);
    } else {
      p_ = p;
    }
  }

  void reset() { bin.reset(); }
  _uint n() { return n_; }
  double p() { return p_; }

  template <class Generator>
  _uint operator()(Generator& g) {
    _uint num_ones = bin(g);
    _uint val = 0;

    for (_uint j = n_ - num_ones; j < n_; ++j) {
      std::uniform_int_distribution<_uint> unif(0, j);
      _uint t = 1 << unif(g);
      if ((val & t) != 0) {
        val = val | (1 << j);
      } else {
        val = val | t;
      }
    }
    return invert ^ val;
  }
};

//implementation found in Watanabe et al.
class PoissonOr {
private:
  _uint n_;
  double p_;
  double lambda;
  std::poisson_distribution<_uint> poiss;
  std::uniform_int_distribution<_uint> unif;
  _uint inversion_layer;
  _uint invert = 0;

public:
  PoissonOr(_uint n, double p) : unif(0, n-1) {
    n_ = n;
    p_ = p;
    if (p > 0.5) {
      p = 1 - p;
      invert = (1 << n) - 1;
    }
    lambda = -(double)n*log(1-p);
    poiss = std::poisson_distribution<_uint>(lambda);
  }

  void reset() { poiss.reset(); }
  _uint n() { return n_; }
  double p() { return p_; }

  template <class Generator>
  _uint operator()(Generator& g) {
    _uint n_strings = poiss(g);
    _uint val = 0;
    for (_uint j = 0; j < n_strings; ++j) {
      val = val | (1 << unif(g));
    }
    return invert ^ val;
  }
};

//binomial shuffler using 2 rv generation
class BinomialShuffleNew {
private:
  _uint n_;
  double p_;
  std::binomial_distribution<_uint> bin;
  std::vector<_uint> choose_vals;
  _uint invert = 0;

  _uint get_bit_stream(_uint n, _uint k, _uint x) {
    _uint ret = 0;
    //correction factor to get in terms of n-1 choose k-1
    _uint chooseval = choose_vals[k]*k/n;
    //note that we only go to i=n-2 to avoid a divide by zero, the i=n-1 case is handled after the loop
    _uint i = 0;
    for (; k > 0 && i < n-1; ++i) {
      //std::cout << "choose: " << choose((n-i)-1, k-1) << ", nval: " << chooseval << std::endl;
      if (x < chooseval ) {
        ret = ret | (0x01 << i);
        k--;
        //special transformation to map chooseval to choose(n-(i+1)-1, (k-1)-1)
        chooseval *= k;
        chooseval /= n-i-1;
      } else {
        x -= chooseval;
        //special transformation to map chooseval to choose(n-(i+1)-1, k-1)
        chooseval *= n-i-k;
        chooseval /= n-i-1;
      }
    }
    if (x < chooseval) {
      ret = ret | 1 << i;
    }
    return ret;
  }

public:
  BinomialShuffleNew(_uint n, double p) : choose_vals(n+1), bin(n, p) {
    n_ = n;
    p_ = p;
    if (p > 0.5) {
      p_ = 1 - p;
      invert = (1 << n) - 1;
      bin = std::binomial_distribution<_uint>(n, p_);
    } else {
      p_ = p;
    }
    
    for (_uint i = 0; 2*i <= n; ++i) {
      choose_vals[i] = choose(n, i);
      choose_vals[n-i] = choose_vals[i];
    }
  }

  void reset() { bin.reset(); }
  _uint n() { return n_; }
  double p() { if (invert) { return 1 - p_; } else { return p_; } }

  template <class Generator>
  _uint operator()(Generator& g) {
    _uint num_ones = bin(g);

    if (num_ones == 0) {
      return invert;
    }
    //minus 1 because we start indexing from 0
    std::uniform_int_distribution<_uint> unif(0, choose_vals[num_ones] - 1);
    _uint bin_val = get_bit_stream(n_, num_ones, unif(g));
    return invert ^ bin_val;
  }
};

//binomial shuffler using 2 rv generation
class BinomialShufflePrecompute {
private:
  double prob_1 = 0.0;
  _uint n_samples = 0;

  _uint n_;
  double p_;
  std::binomial_distribution<_uint> bin;
  std::vector<_uint> choose_vals;
  const _uint group_size = 32;
  _uint group_n = group_size;
  _uint n_groups = 1;
  //const _uint int byte_size = 8;
  _uint precompute_k;
  std::vector<std::vector<_uint>> precompute_strings;
  _uint invert = 0;
  _uint mask = 0;

  _uint get_bit_stream(_uint n, _uint k, _uint x, _uint precompute_limit=0) {
    _uint ret = 0;
    //correction factor to get in terms of n-1 choose k-1
    _uint chooseval = choose_vals[k]*k/n;
    //note that we only go to i=n-2 to avoid a divide by zero, the i=n-1 case is handled after the loop
    _uint i = 0;
    for (; k > precompute_limit && i < n-1; ++i) {
      //std::cout << "choose: " << choose((n-i)-1, k-1) << ", nval: " << chooseval << std::endl;
      if (x < chooseval ) {
        ret = ret | (0x01 << i);
        k--;
        //special transformation to map chooseval to choose(n-(i+1)-1, (k-1)-1)
        chooseval *= k;
        chooseval /= n-i-1;
      } else {
        x -= chooseval;
        //special transformation to map chooseval to choose(n-(i+1)-1, k-1)
        chooseval *= n-i-k;
        chooseval /= n-i-1;
      }
    }
    if (x < chooseval) {
      ret = ret | 1 << i;
      k--;
    }
    if (k > 0 && k <= precompute_limit) {
      ret |= precompute_strings[k][choose_vals[k]-x-1];
    }
    return ret;
  }

public:
  BinomialShufflePrecompute(_uint n, double p, double p_precompute_k=4) : choose_vals(n+1), bin(n, p), precompute_strings(p_precompute_k) {
    n_ = n;
    p_ = p;
    //n_bytes = (n + (byte_size-1))/n;
    //for large bitstrings, it is helpful to chunk results, this mask contains n_ one bits and is safe for n_=<word size>
    mask = ( ((_uint)1 << (n_-1)) - 1 ) | ( (_uint)1 << (n_-1) );
    n_groups = (n + (group_size-1))/group_size;
    
    if (n > group_size) {
      n = group_size;
      bin = std::binomial_distribution<_uint>(n, p_);
    }
    group_n = n;
    precompute_k = p_precompute_k;
    
    if (p > 0.5) {
      p_ = 1 - p;
      invert = mask;
      bin = std::binomial_distribution<_uint>(n, p_);
    } else {
      p_ = p;
    }
    
    for (_uint i = 0; 2*i <= group_n; ++i) {
      choose_vals[i] = choose(group_n, i);
      choose_vals[n-i] = choose_vals[i];
    }

    for (_uint k = 0; k < precompute_k; ++k) {
      _uint chooseval = choose(group_n, k);
      precompute_strings[k].resize(chooseval);
      for (_uint j = 0; j < chooseval; ++j) {
        precompute_strings[k][j] = get_bit_stream(group_n, k, j);
      }
    }
  }

  void reset() { bin.reset(); }
  _uint n() { return n_; }
  double p() { if (invert) { return 1 - p_; } else { return p_; } }

  template <class Generator>
  _uint operator()(Generator& g) {
    _uint result = 0;
    for (_uint i = 0; i < n_groups; ++i) {
      _uint num_ones = bin(g);

      if (num_ones != 0) {
        //minus 1 because we start indexing from 0
        std::uniform_int_distribution<_uint> unif(0, choose_vals[num_ones] - 1);
        _uint j = unif(g);
        if (num_ones < precompute_k) {
          result |= precompute_strings[num_ones][j] << (i*group_size);
        } else {
          result |= get_bit_stream(group_n, num_ones, j, precompute_k-1) << (i*group_size);
        }
      }
    }
    if ( 1 & (invert ^ result) ) {
      prob_1 = (prob_1*(double)n_samples + 1.0)/(n_samples+1);
    } else {
      prob_1 = prob_1*(double)n_samples / (n_samples+1);
    }
    ++n_samples;
    if ( (1 << (n_-1)) & (invert ^ result) ) {
      prob_1 = (prob_1*(double)n_samples + 1.0)/(n_samples+1);
    } else {
      prob_1 = prob_1*(double)n_samples / (n_samples+1);
    }
    ++n_samples;
    if (n_samples % 10000 == 0) {
      std::cout << "prob_" << n_samples << ":" << prob_1 << std::endl;
    }
    return mask & (invert ^ result);
  }
};

#endif //BITSTRINGS_H
