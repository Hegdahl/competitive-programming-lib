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
  enum Layer {
    ST, L, M1, M2, R
  };
  
  int S = 0;
  int T = 1;
  Flow<int> flow(2, n, n, n, n);
 
  auto add_edge = [&](int x, int y) {
    int i = where[x];
    int j = where[y];
    if (i == -1 || j == -1) return;
    flow(L, i)(M1, j) = 1;
    flow(M2, i)(R, j) = 1;
  };
 
  for (int i = 0; i < n; ++i) {
    flow(S)(L, i) = 1;
    flow(M1, i)(M2, i) = 1;
    flow(R, i)(T) = 1;
 
    if (a[i] != 1)
      add_edge(a[i], 1);
    for (int d = 2; d*d <= a[i]; ++d) {
      if (a[i] % d) continue;
 
      add_edge(a[i], d);
      if (d*d != a[i])
        add_edge(a[i], a[i]/d);
    }
  }
 
  std::cout << flow(S)(T).max_flow() << '\n';
 
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
