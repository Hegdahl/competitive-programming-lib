#pragma once
#include <number_theory/factor_sieve.hpp>

template <int N>
struct Mobius {
  int mu[N + 1];

  constexpr Mobius(const FactorSieve<N> &sieve) {
    mu[0] = 0;
    mu[1] = 1;
    for (int i = 2; i <= N; ++i) {
      if (sieve.smallest_factor[i] ==
          sieve.smallest_factor[i / sieve.smallest_factor[i]])
        mu[i] = 0;
      else
        mu[i] = -mu[i / sieve.smallest_factor[i]];
    }
  }

  int operator[](size_t i) const { return mu[i]; }
};