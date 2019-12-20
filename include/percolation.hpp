#ifndef PERCOLATION_H
#define PERCOLATION_H

#include "bitstrings.hpp"

template <typename Sampler>
class DirectedPercolator1 {
private:
  _uint string_size;
  _uint last_bit;
  _uint t_max;
  _uint t = 0;
  _uint steady_state_t;
  double p;
//  std::vector<std::vector<_uint>> bitstrings;//bitstrings[time][spatial location]
  std::vector<_uint> bitstring;

  _uint one_state(_uint t = 0) {
    if (t == 0) { t = string_size; }
    if (t >= string_size) { return (last_bit-1) | last_bit; }
    return (1 << t) - 1;
  }
  _uint n_samples = 0;
  double test_p0, test_pl, test_pm;

public:
  DirectedPercolator1(double p_p, _uint pt_max, bool relax, _uint p_string_size) : bitstring(pt_max + 1, 0) {
    t_max = pt_max;
    steady_state_t = t_max;
    string_size = p_string_size;
    last_bit = 1 << (string_size - 1);
    /*bitstrings.resize(t_max);
    _uint n_strings = 0;
    for (_uint i = 0; i < t_max; ++i) {
      if (i % string_size == 0) {
        ++n_strings;
      }
      bitstrings[i].resize(n_strings);
    }
    bitstrings[0][0] = 1;*/
    if (relax) {
      for (_uint i = 0 ; i < bitstring.size(); ++i) {
        //set everything to be the all ones state
        bitstring[i] = last_bit | (last_bit-1);
      }
    } else {
      bitstring[0] = 1;
    }
    test_p0 = 0;
    test_pl = 0;
    test_pm = 0;
  }

  template <typename Generator>
  void update(Sampler& s, Generator& g, bool relax) {
    _uint mask = last_bit | (last_bit - 1);
    ++t;
    std::vector<_uint> old_bitstring = bitstring;
    _uint n_strings = (t + string_size) / string_size;
    //std::vector<_uint> l_bonds(n_strings);
    //std::vector<_uint> r_bonds(n_strings);
    _uint l_bonds;
    _uint r_bonds;
    if (t < steady_state_t && t < t_max) {
      bool all_zeros = true;

      for (_uint i = 0; i < n_strings; ++i) {
        //handle carryover
        if (i > 0 && (old_bitstring[i-1] & last_bit) > 0) {
          if ((l_bonds & last_bit) > 0) {
            bitstring[i] = 1;
          } else {
            bitstring[i] = 0;
          }
        } else {
          bitstring[i] = 0;
        }
        //perform updates
        l_bonds = s(g);
        r_bonds = s(g);

        bitstring[i] |= (old_bitstring[i] & r_bonds);
        bitstring[i] |= ( (old_bitstring[i] & l_bonds) << 1 ) & mask;
        if (bitstring[i] != 0) { all_zeros = false; }
      }
      if (all_zeros) {
        steady_state_t = t;
      }
    }
  }

  template <typename Generator>
  void update_debug(Sampler& s, Generator& g, bool relax) {
    _uint mask = last_bit | (last_bit - 1);
    if (t == 0) {
      if (relax) {
        for (_uint i = 0; i < t_max; ++i) {
          std::cout << "1 ";
        }
        std::cout << "1\n";
      } else {
        for (_uint i = 0; i < t_max - (t+1); ++i) {
          std::cout << " ";
        }
        std::cout << "1\n";
      }
    }
    ++t;
    std::vector<_uint> old_bitstring = bitstring;
    _uint n_strings = (t + string_size) / string_size;
    std::vector<_uint> l_bonds(n_strings);
    std::vector<_uint> r_bonds(n_strings);
    if (t < steady_state_t && t < t_max) {
      bool all_zeros = true;

      for (_uint i = 0; i < n_strings; ++i) {
        //handle carryover
        if (i > 0 && (old_bitstring[i-1] & last_bit) > 0) {
          if ((l_bonds[i-1] & last_bit) > 0) {
            bitstring[i] = 1;
          } else {
            bitstring[i] = 0;
          }
        } else {
          bitstring[i] = 0;
        }
        //perform updates
        //_uint l_bonds = s(g);
              //_uint r_bonds = s(g);
        l_bonds[i] = s(g);
        r_bonds[i] = s(g);
        _uint samples0 = (l_bonds[i] & 1) + (r_bonds[i] & 1);
        _uint samplesm = ((l_bonds[i] >> (string_size/2)) & 1) + ((r_bonds[i] >> (string_size/2)) & 1);
        _uint samplesl = ((l_bonds[i] >> (string_size-1)) & 1) + ((r_bonds[i] >> (string_size-1)) & 1);
        test_p0 = ( test_p0*n_samples + (double)(samples0) )/(n_samples+2);
        test_pm = ( test_pm*n_samples + (double)(samplesm) )/(n_samples+2);
        test_pl = ( test_pl*n_samples + (double)(samplesl) )/(n_samples+2);
        n_samples += 2;

        bitstring[i] |= (old_bitstring[i] & r_bonds[i]);
        bitstring[i] |= ( (old_bitstring[i] & l_bonds[i]) << 1 ) & mask;
        if (bitstring[i] != 0) { all_zeros = false; }
        std::cout << bitstring[i] << " ";
      }
      std::cout << "\n";
      std::cout << "p0: " << test_p0 << ", pm" << test_pm << ", pl" << test_pl << std::endl;
      if (all_zeros) {
        steady_state_t = t;
      }
      _uint old_t = t;
      if (relax) {
        t = t_max;
      } else {
        for (_uint i = 1; i < t_max - (t+1); ++i) {
          std::cout << " ";
        }
      }

      _uint status = 0;
      for (_uint i = 0; i <= t; ++i) {
        _uint ind = (t-i) / string_size;
        _uint bit = (t-i) % string_size;
        status = ( 2 & (l_bonds[ind] >> (bit-1)) ) | ( 1 & (r_bonds[ind] >> bit) ) ;
        if ( ((old_bitstring[ind] >> bit) & 1) != 0 && i < t) {
          if (status == 0) {
            std::cout << " ";
          } else if (status == 1) {
            std::cout << "\\";
          } else if (status == 2) {
            std::cout << "/";
          } else {
            std::cout << "^";
          }
        } else {
          std::cout << " ";
        }
        if (i < t) {
          std::cout << " ";
        }
      }
      if ((old_bitstring[0] & 1) != 0) {
        status = (r_bonds[0] & 1) | ((l_bonds[0] & 1)<< 1);
        if (status == 0) {
          std::cout << " ";
        } else if (status == 1) {
          std::cout << "\\";
        } else if (status == 2) {
          std::cout << "/";
        } else {
          std::cout << "^";
        }
      }
      std::cout << "\n";
      for (_uint i = 1; i < t_max - t; ++i) {
        std::cout << " ";
      }
      for (_uint i = 0; i <= t; ++i) {
        _uint ind = (t-i) / string_size;
        _uint bit = (t-i) % string_size;
        std::cout << ( (bitstring[ind] >> bit) & 1 ) << " ";
        _uint status = ( 2 & (l_bonds[ind] >> (bit-1)) ) | ( 1 & (r_bonds[ind] >> bit) );
      }
      std::cout << "\n";
      t = old_t;
    }
  }

  _uint get_n_t() {
    if (t > steady_state_t) { return 0; }

    _uint n_strings = (t + string_size) / string_size;
    _uint count = 0;
    for (_uint i = 0; i < n_strings; ++i) {
      _uint tmp = bitstring[i];
      while (tmp > 0) {
	count += tmp & 1;
	tmp = tmp >> 1;
      }
    }
    return count;
  }

  double get_rho_t() {
    return (double)(get_n_t())/(t+1);
  }

  std::vector<bool> get_state() {
    std::vector<bool> result;
    result.reserve(t+1);
    for (_uint i = 0; i < t+1; ++i) {
      _uint ind = i / string_size;
      _uint bit = i % string_size;
      result.push_back((bitstring[ind] >> bit) & 1);
    }
    return result;
  }
  
  bool alive() {
    if (t > steady_state_t) {
      return false;
    }
    return true;
  }
};

template <typename Sampler>
class PercolationTracker {
private:
  std::vector<double> avg_occupation_n;
  std::vector<double> avg_rho;
  std::vector<_uint> n_survivors;
  Sampler s;
  std::vector<DirectedPercolator1<Sampler>> percs;
  _uint string_size;
  bool relax = false;
  _uint t_max;

public:
  PercolationTracker(_uint n_samples, double p, bool p_relax=false, _uint p_string_size=32, _uint pt_max=1000) :
  s(p_string_size, p),
  avg_occupation_n(1),
  avg_rho(1),
  n_survivors(1),
  relax(p_relax),
  t_max(pt_max) {
    /*if (!s.bijectivity_test()) {
      std::cout << "ERROR\n";
    }*/
    string_size = p_string_size;

    percs.reserve(n_samples);
    for (_uint i = 0; i < n_samples; ++i) {
      percs.emplace_back(p, pt_max, relax, string_size);
    }
    avg_occupation_n[0] = 1;
    avg_rho[0] = 1;
    n_survivors[0] = n_samples;
  }

  template <typename Generator>
  void update(Generator& g) {
    _uint t = avg_occupation_n.size();
    avg_occupation_n.push_back(0);
    avg_rho.push_back(0);
    n_survivors.push_back(0);
    for (_uint i = 0; i < percs.size(); ++i) {
      percs[i].update(s, g, relax);
      //information for debugging
      /*std::vector<bool> status = percs[i].get_state();
      for(_uint j = 0; j < status.size(); ++j) {
        std::cout << status[j];
      }
      std::cout << std::endl;*/
      avg_occupation_n[t] += (double)(percs[i].get_n_t()) / percs.size();
      if (relax) {
        avg_rho[t] += (double)(percs[i].get_n_t()) / (percs.size()*t_max);
      } else {
        avg_rho[t] += percs[i].get_rho_t() / percs.size();
      }
      if (percs[i].alive()) {
        ++n_survivors[t];
      }
    }
  }

  std::vector<double> get_time_arr() {
    std::vector<double> time(avg_rho.size());
    for (_uint i = 0; i < time.size(); ++i) {
      time[i] = (double)i;
    }
    return time;
  }
  std::vector<double> get_occupation_n_arr() { return avg_occupation_n; }
  std::vector<double> get_rho_arr() { return avg_rho; }
  std::vector<_uint> get_survivors_arr() { return n_survivors; }

  double get_n(_uint t) { return avg_occupation_n[t]; }
  double get_rho(_uint t) { return avg_rho[t]; }
  double get_survivors(_uint t) { return n_survivors[t]; }
};

#endif //PERCOLATION_H
