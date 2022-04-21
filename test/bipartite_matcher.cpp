/**
 * problem: Matching on Bipartite Graph
 * url:     https://judge.yosupo.jp/problem/bipartitematching
 */
#include <graph/bipartite_matcher.hpp>

#include <iostream>

int main() {
  std::cin.tie(0)->sync_with_stdio(0);

  int L, R, M;
  std::cin >> L >> R >> M;
  BipartiteMatcher bm(L, R);

  while (M--) {
    int i, j;
    std::cin >> i >> j;
    bm.add_edge(i, j);
  }

  std::cout << bm.compute_matching_size() << '\n';
  return 0;

  auto matching = bm.get_matching();
  std::cout << matching.size() << '\n';
  for (auto [i, j] : matching)
    std::cout << i << ' ' << j << '\n';
}
