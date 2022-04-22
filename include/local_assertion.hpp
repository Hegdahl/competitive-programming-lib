#pragma once

#include <iostream>

#ifdef ENABLE_DEBUG
#define LOCAL_ASSERT(condition, message)            \
  if (!(condition)) {                               \
    std::cerr << "Assertion (" << #condition        \
              << ") failed\non line " << __LINE__   \
              << " in " << __FILE__                 \
              << '#' << __FUNCTION__ << "()\n";     \
    std::cerr << message << '\n';                   \
    abort();                                        \
  }
#else
#define LOCAL_ASSERT(...) do {} while (0);
#endif

#define LOCAL_ASSERT_IN_RANGE(variable, low, high)  \
  LOCAL_ASSERT(low <= variable && variable <= high, \
               #variable << '=' << variable         \
               << " is out of the range "           \
               << '[' << low << ", " << high << ']')
