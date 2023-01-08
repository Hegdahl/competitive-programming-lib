#pragma once

#include <ints.hpp>
#include <type_traits>
#include <vector>

/**
 * @brief Maintain prefix sums of conceptual array A.
 *
 * Construction is O(n), and all operations are O(log n).
 *
 * @tparam T type of values to hold. Should be an abelian group
 *         using operators + and - where default construction
 *         results in the identity element.
 */
template <class T, class index_t = u32,
          i32 max_log2n = 8 * sizeof(index_t) - bool(std::is_signed_v<index_t>),
          class Container = std::vector<T>, bool init_loop = true>
class FenwickTree {
  Container bit;

 public:
  /**
   * @brief Construct a new Fenwick Tree object
   * as if it was A was a std::vector<T>.
   *
   * @param args arguments to forward to the container constructor
   */
  template <class... Args>
  FenwickTree(Args &&...args) : bit(std::forward<Args>(args)...) {
    if constexpr (init_loop) {
      for (index_t i = 0; i != size(); ++i) {
        index_t j = i | (i + 1);
        if (j < size()) bit[j] += bit[i];
      }
    }
  }

  /**
   * @brief A[idx] += x.
   *
   * @param idx index of value to add to
   * @param x   what to add to the value
   */
  void update(index_t idx, const T &x) {
    while (idx < size()) {
      bit[idx] += x;
      idx |= idx + 1;
    }
  }

  /**
   * @brief Find Σ(A[:end]).
   *
   * @param end beyond last idx that is part of the sum
   * @return T  sum of values with idx smaller than end
   */
  T operator()(index_t end) const {
    T x{};
    while (end) {
      x += bit[end - 1];
      end &= end - 1;
    }
    return x;
  }

  /**
   * @brief Find largest idx such that Σ(A[:idx]) <= k.
   *
   * Assumes prefix sums are monotonic.
   *
   * @param k max allowed prefix sum
   * @return (idx, k - Σ(A[:idx]))
   */
  std::pair<index_t, T> max_leq(T k) const {
    index_t idx = -index_t(1);
    for (i32 d = max_log2n - 1; d >= 0; --d) {
      index_t right_idx = idx + (index_t(1) << d);
      if (right_idx < size() && bit[right_idx] <= k) {
        idx = right_idx;
        k -= bit[idx];
      }
    }
    return {idx + 1, k};
  }

  index_t size() const { return (index_t)bit.size(); }
};
