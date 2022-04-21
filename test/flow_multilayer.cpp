/**
 * problem: F. Making It Bipartite
 * url:     https://codeforces.com/contest/1630/problem/F
 */
#include <graph/flow.hpp>
#include <iostream>

std::vector<int> where((int)5e4+1, -1);
 
void tcase() {
  int n;
  std::cin >> n;
  std::vector<int> a(n);
 
  for (int &x : a)
    std::cin >> x;
 
  for (int i = 0; i < n; ++i)
    where[a[i]] = i;
 
  // s l m1 m2 r t
  
  Flow<int> flow(1, n, n, n, n, 1);
 
  auto add_edge = [&](int x, int y) {
    int i = where[x];
    int j = where[y];
    if (i == -1 || j == -1) return;
    flow(1, i)(2, j) = 1;
    flow(3, i)(4, j) = 1;
  };
 
  for (int i = 0; i < n; ++i) {
    flow(0, 0)(1, i) = 1;
    flow(2, i)(3, i) = 1;
    flow(4, i)(5, 0) = 1;
 
    if (a[i] != 1)
      add_edge(a[i], 1);
    for (int d = 2; d*d <= a[i]; ++d) {
      if (a[i] % d) continue;
 
      add_edge(a[i], d);
      if (d*d != a[i])
        add_edge(a[i], a[i]/d);
    }
  }
 
  std::cout << flow(0, 0)(5, 0).max_flow() << '\n';
 
  for (int x : a)
    where[x] = -1;
}
 
int main() {
  std::cin.tie(0)->sync_with_stdio(0);
 
  int T;
  std::cin >> T;
  while (T--)
    tcase();
}
