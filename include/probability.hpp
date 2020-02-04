#include <math.h>
#include <random>
#include <vector>
#include <iostream>

typedef size_t _uint;

struct Range {
  double min;
  double max;
  Range(double p_min=0, double p_max=0) : min(p_min), max(p_max) {}
  double avg() { return (max+min) / 2; }
  double gap() { return (max-min); }
};

double choose_approx( _uint n, _uint k ) {
  if (n > 64) {
    double ratio = (double)k/n;
    double term1 = pow(1/ratio, k);
    double term2 = pow(1/(1 - ratio), n-k);
    double term3 = sqrt(6.283186*ratio*(n-k));
    term1 = term1 / term3;
    term2 = term2;
    return term1*term2;
  } else {
    if (k > n) return 0;
    if (k * 2 > n) k = n-k;
    if (k == 0) return 1;

    _uint result = n;
    for( _uint i = 2; i <= k; ++i ) {
      result *= (n-i+1);
      result /= i;
    }
    return (double)result;
  }
}

//stolen from https://math.stackexchange.com/questions/1313831/inverse-complementary-error-function-values-near-0
double erf_inv(double x) {
  double p_x = x - 4397*M_PI*pow(x, 3)/17352 + 111547*M_PI*M_PI*pow(x, 5)/14575680;
  double q_x = 1 - 5843*M_PI*pow(x, 2)/17352 + 20533*M_PI*M_PI*pow(x, 4)/971712;
  return sqrt(M_PI)*p_x / q_x;
}

Range quadratic(double a, double b, double c) {
  Range result;
  result.min = (-b - sqrt(b*b - 4*a*c)) / (2*a);
  result.max = (-b + sqrt(b*b - 4*a*c)) / (2*a);
  return result;
}

template<typename Sampler>
class BitsData {
private:
  std::vector<unsigned> data;
  Sampler s;
  _uint n_trials = 0;

  //n = number of observed ones, m = number of trials
  Range find_uncertainty(unsigned n, unsigned m, unsigned iterations=5, double confidence = 0.95) {
    //checks to avoid division by zero
    if (n == 0) {
      return Range(0, 1 - pow(confidence, 1.0/m));
    }
    if (n >= m) {
      return Range(pow(confidence, 1.0/m), 1);
    }

    double chi = 2*(double)m*pow(erf_inv(2*confidence - 1), 2);
    double u_obs = (double)n;

    Range result = quadratic(m*m + chi, -2*(double)m*u_obs-chi, u_obs*u_obs);
    if (result.max > 1) { result.max = 1; }
    if (result.min > 1) { result.min = 1; }
    if (result.max < 0) { result.max = 0; }
    if (result.min < 0) { result.min = 0; }
    return result;
  }
public:
  BitsData(_uint p_len, double p_p) : s(p_len, p_p), data(p_len, 0) {}

  void update(std::mt19937& generator) {
    _uint sample = s(generator);
    for (_uint j = 0; j < data.size(); ++j) {
      if ( ((sample >> j) & 1) != 0 ) {
        data[j] += 1;
      }
    }
    ++n_trials;
  }

  void print_statistics(double p_theory) {
    for (_uint j = 0; j < data.size(); ++j) {
      Range r = find_uncertainty(data[j], n_trials);
      std::cout << "\tp_" << j << "=" << r.avg() << "\u00B1" << r.gap()/2;
      if (r.max < p_theory || r.min > p_theory) {
        std::cout << "!";
      }
      std::cout << "\n";
    }
  }
};
