#pragma once

#include <array>
#include <new>

/**
 * @brief Constant optimized non-lazy segment tree
 *
 * Mostly taken from
 * https://en.algorithmica.org/hpc/data-structures/segment-trees/
 */
template <class Group, int N>
struct SIMDSegTree {
  static_assert(Group::commute);
  using T = typename Group::value_type;

  static constexpr int cache_size = 64; // hope it's 64 for clang support... std::hardware_constructive_interference_size;
  static constexpr int reg_size = 32;  // assuming avx2
  static constexpr int reg_count = cache_size / reg_size;
  static constexpr int branching = cache_size / sizeof(T);
  static constexpr int branching_bits = std::bit_width(unsigned(branching)) - 1;

  // MUST BE TYPEDEF, NOT USING
  typedef T __attribute__((vector_size(reg_size))) vec;

  struct Precalc {
    alignas(cache_size) T mask[branching][branching];

    constexpr Precalc() : mask{} {
      for (int i = 0; i < branching; ++i)
        for (int j = 0; j < branching; ++j) mask[i][j] = (i < j ? -1 : 0);
    }
  };

  static constexpr Precalc precalc{};

  static constexpr int calc_height(int n) {
    return (n <= branching ? 1 : calc_height(n / branching) + 1);
  }

  static constexpr int calc_offset(int h) {
    int s = 0, n = N;
    while (h--) {
      s += (n + branching - 1) / branching * branching;
      n /= branching;
    }
    return s;
  }

  static constexpr int round(int k) {
    return k & ~(branching - 1);  // = k / branching * branching
  }

  static constexpr int H = calc_height(N);

  template<std::size_t ... indices>
  static constexpr std::array<int, H+1> calc_offsets(std::index_sequence<indices...>) {
    return {calc_offset(indices)...};
  }

  static constexpr auto offset = calc_offsets(std::make_index_sequence<H+1>());

  alignas(cache_size) T values[offset[H]]{};

  void add(int k, T x) {
    vec v = x + vec{};
    for (int h = 0; h != H; ++h) {
      auto a = (vec *)&values[offset[h] + round(k)];
      const auto m = (const vec *)&precalc.mask[k % branching];
      for (int i = 0; i != reg_count; ++i) a[i] = Group::op(a[i], v & m[i]);
      k >>= branching_bits;
    }
  }

  T sum(int i) {
    T s = 0;
    for (int h = 0; h != H; ++h) {
      s = Group::op(s, values[offset[h] + i]);
      i >>= branching_bits;
    }
    return s;
  }

  T sum(int i, int j) {
    T s = 0;
    for (int h = 0; h != H; ++h) {
      s = Group::op(s, Group::inverse(values[offset[h] + i]));
      s = Group::op(s, values[offset[h] + j]);
      i >>= branching_bits;
      j >>= branching_bits;
    }
    return s;
  }
};

template<class T>
struct Add {
  using value_type = T;
  static constexpr bool commute = true;
  static constexpr auto op(auto &&a, auto&&b) { return a + b; }
  static constexpr auto inverse(auto &&a) { return -a; }
  static constexpr auto repeat(auto &&a, long long n) { return a * n; }
  static constexpr T id() { return 0; }
};

template<class T>
struct Multiply {
  using value_type = T;
  static constexpr bool commute = true;
  static constexpr auto op(auto &&a, auto&&b) { return a * b; }
  static constexpr auto inverse(auto &&a) { return 1 / a; }
  static constexpr T id() { return 1; }
};

template<class T>
struct Xor {
  using value_type = T;
  static constexpr bool commute = true;
  static constexpr auto op(auto &&a, auto&&b) { return a ^ b; }
  static constexpr auto inverse(auto &&a) { return a; }
  static constexpr auto repeat(auto &&a, long long n) { return a * (n & 1); }
  static constexpr T id() { return 0; }
};

