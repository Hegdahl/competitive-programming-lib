#pragma once

#define _ __unused __attribute__((unused))
#include <impyster/containers.hpp>
#include <impyster/enumerate.hpp>
#include <impyster/map.hpp>
#include <impyster/print.hpp>
#include <impyster/range.hpp>
#include <impyster/zip.hpp>
#include <iostream>

namespace impyster {

struct {
  template<class T>
  operator T() {
    T x;
    std::cin >> x;
    return x;
  }
} in;

} // namespace impyster
