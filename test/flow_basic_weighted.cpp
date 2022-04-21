/**
 * problem: Download Speed
 * url:     https://cses.fi/problemset/task/1694
 */

#include <graph/flow.hpp>

#include <iostream>

int main() {
  std::cin.tie(0)->sync_with_stdio(0);

  int n, m;
  std::cin >> n >> m;
  Flow f(n);
  while (m--) {
    int i, j, w;
    std::cin >> i >> j >> w;
    --i, --j;
    f(i)(j) = w;
  }

  std::cout << f(0)(n-1).max_flow() << '\n';
}
