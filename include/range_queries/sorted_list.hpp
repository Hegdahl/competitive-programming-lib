#pragma once

#include <ints.hpp>
#include <local_assertion.hpp>
#include <range_queries/fenwick_tree.hpp>

#include <algorithm>
#include <array>
#include <bitset>
#include <type_traits>

template <class T, class Cmp = std::less<T>, u16 max_block_size = 256>
class SortedList {

  struct Block {
    using index_t = std::conditional_t<max_block_size <= 256, u8, u16>;
    static constexpr u32 max_jump = []() {
      u32 res = 1;
      while ((res << 1) <= max_block_size) res <<= 1;
      return res;
    }();
    std::array<T, max_block_size> values;
    std::array<index_t, max_block_size> indices;
    std::bitset<max_block_size> alive;
    u32 size = 0;

    index_t next_available() const {
      LOCAL_ASSERT(alive.count() < max_block_size,
                   "Tried find available in full block");
      static index_t circular_pointer = 0;
      while (alive[circular_pointer])
        circular_pointer = (circular_pointer + 1) % max_block_size;
      return circular_pointer;
    }

    u32 lower_bound(const T &x, const Cmp &cmp) const {
      u32 p = 0;
      for (u32 jump = max_jump; jump; jump >>= 1)
        if (p + jump <= size && cmp(values[indices[p + jump - 1]], x))
          p += jump;
      return p;
    }

    u32 upper_bound(const T &x, const Cmp &cmp) const {
      u32 p = 0;
      for (u32 jump = max_jump; jump; jump >>= 1)
        if (p + jump <= size && !cmp(x, values[indices[p + jump - 1]]))
          p += jump;
      return p;
    }

    template <class U>
    void insert(U &&x, const Cmp &cmp) {
      LOCAL_ASSERT(size < max_block_size, "Tried to insert in full block.");

      u32 p = upper_bound(x, cmp);

      for (u32 q = size; q > p; --q) indices[q] = indices[q - 1];
      indices[p] = next_available();
      alive[indices[p]] = 1;
      values[indices[p]] = std::forward<U>(x);
      ++size;
    }

    void erase(const T &x, const Cmp &cmp) {
      LOCAL_ASSERT(0 < size, "Tried to erase from empty block.");

      u32 p = upper_bound(x, cmp) - 1;
      LOCAL_ASSERT(values[indices[p]] == x,
                   "Tried to erase element that is not in the block.");

      --size;
      alive[indices[p]] = 0;
      for (u32 q = p; q < size; ++q) indices[q] = indices[q + 1];
    }

    T erase_at(u32 p) {
      T res = values[indices[p]];
      --size;
      alive[indices[p]] = 0;
      for (u32 q = p; q < size; ++q) indices[q] = indices[q + 1];
      return res;
    }

    const T &operator[](u32 p) const { return values[indices[p]]; }

    void split(Block &right) {
      LOCAL_ASSERT(size == max_block_size,
                   "Tried to split block without being full.");

      size = (max_block_size + 1) / 2;

      right.size = max_block_size - size;
      for (u32 p = size; p != max_block_size; ++p) {
        right.values[p - size] = std::move(values[indices[p]]);
        right.indices[p - size] = p - size;
        right.alive[p - size] = 1;
        alive[indices[p]] = 0;
      }
    }
  };

  Cmp cmp_;
  size_t size_;
  std::vector<Block> block_;
  std::vector<u32> block_indices_;
  FenwickTree<u32> psum;
  std::vector<T> separators;

 public:
  SortedList(const Cmp &cmp = Cmp())
      : cmp_(cmp),
        size_(0),
        block_(1),
        block_indices_(1),
        psum(1),
        separators() {}

  template <class U>
  void insert(U &&x) {
    int i = block_of(x);

    if (block(i).size == max_block_size) {
      split(i);
      if (x >= separators[i]) ++i;
    }

    block(i).insert(std::forward<U>(x), cmp_);
    psum.update(i, 1);
    ++size_;
  }

  void erase(const T &x) {
    int i = block_of(x);
    if (block(i).size == 0 ||
        (i && separators[i - 1] == x && block(i)[0] != x)) {
      pop(psum(i) - 1);
    } else {
      block(i).erase(x, cmp_);
      psum.update(i, -1);
      --size_;
    }
  }

  T pop(i32 k = 0) {
    if (k < 0) k += size();
    auto [i, j] = psum.max_leq(k);
    T removed = block(i).erase_at(j);
    psum.update(i, -1);
    --size_;
    return removed;
  }

  const T &operator[](i32 k) const {
    if (k < 0) k += size();
    auto [i, j] = psum.max_leq(k);
    return block(i)[j];
  }

  size_t size() const { return size_; }

  T next(const T &x) const {
    u32 i = block_of(x);
    u32 j = block(i).upper_bound(x, cmp_);
    if (j < block(i).size)
      return block(i)[j];
    else
      return (*this)[psum(i+1)];
  }

  T prev(const T &x) const {
    u32 i = earliest_block_of(x);
    u32 j = block(i).lower_bound(x, cmp_);

    if (j)
      return block(i)[j - 1];
    else
      return (*this)[psum(i) - 1];
  }

  void dbg() {
    /*
    std::cerr << "| ";
    for (u32 bi : block_indices_) {
      const auto &block = block_[bi];
      for (u32 i = 0; i != block.size; ++i)
        std::cerr << block[i] << ' ';
      std::cerr << "| ";
    }
    std::cerr << '\n'; // */
  }

 private:
  Block &block(u32 i) { return block_[block_indices_[i]]; }
  const Block &block(u32 i) const { return block_[block_indices_[i]]; }

  u32 block_of(const T &x) const {
    u32 i = 0;
    u32 mx = separators.size();
    for (u32 jump = 1u << 28; jump; jump >>= 1)
      if (i + jump <= mx && !cmp_(x, separators[i + jump - 1])) i += jump;
    return i;
  }

  u32 earliest_block_of(const T &x) const {
    u32 i = 0;
    u32 mx = separators.size();
    for (u32 jump = 1u << 28; jump; jump >>= 1)
      if (i + jump <= mx && cmp_(separators[i + jump - 1], x)) i += jump;
    if (i < mx) {
      //std::cerr << "x = " << x << '\n';
      //std::cerr << "separators[i] = " << separators[i] << '\n';
    }
    if (i < mx && separators[i] == x) ++i;
    return i;
  }

  void split(u32 i) {
    block_indices_.insert(block_indices_.begin() + i + 1, block_.size());
    auto &new_block = block_.emplace_back();
    block(i).split(new_block);
    separators.insert(separators.begin() + i, block_.back()[0]);

    psum = FenwickTree<u32>([&]() {
      std::vector<u32> sizes(block_.size());
      for (u32 p = 0; p < sizes.size(); ++p) sizes[p] = block(p).size;
      return sizes;
    }());
  }
};
