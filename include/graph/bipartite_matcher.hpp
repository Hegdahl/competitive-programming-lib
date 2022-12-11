#include <array>
#include <graph/flow.hpp>
#include <local_assertion.hpp>
#include <vector>

/**
 * @brief computes the following :
 *    maximum matching on bipartite graph
 *    minimum vertex cover of bipartite graph
 */
struct BipartiteMatcher {
  Flow<int, int, true> internal_flow;
  int n_l, n_r;

  BipartiteMatcher(int n_l, int n_r)
      : n_l(n_l), n_r(n_r), internal_flow(1, n_l, n_r, 1) {
    for (int i = 0; i < n_l; i++) internal_flow(0, 0)(1, i) = 1;
    for (int i = 0; i < n_r; i++) internal_flow(2, i)(3, 0) = 1;
  }

  BipartiteMatcher(int n) : BipartiteMatcher(n, n) {}

  void add_edge(int u, int v) { internal_flow(1, u)(2, v) = 1; }

  int compute_matching_size() { return internal_flow(0, 0)(3, 0).max_flow(); }

  std::vector<std::array<int, 2>> get_matching() {
    int matching_size = compute_matching_size();

    std::vector<std::array<int, 2>> matching(matching_size);
    auto it = matching.begin();

    for (const auto &pack : internal_flow.get_edge_group(1, 2)) {
      auto &u = std::get<0>(pack);
      auto &v = std::get<1>(pack);
      auto &flow = std::get<2>(pack);
      if (flow > 0) *(it++) = std::array<int, 2>({u, v});
    }

    for (auto &e : matching)
      for (auto &v : e) v = internal_flow.node_name_translator.from_id(v)[1];

    return matching;
  }

  std::vector<int> min_vertex_cover() {
    auto M = get_matching();
    std::vector<int> cover(M.size());
    std::vector<uint8_t> vis(n_l + n_r);
    for (auto &e : M) vis[e[0]] = 1;
    std::vector<std::vector<int>> g(n_l + n_r);
    for (auto &pack : internal_flow.get_edge_group(1, 2)) {
      auto &u = std::get<0>(pack);
      auto &v = std::get<1>(pack);
      auto &flow = std::get<2>(pack);
      if (flow)
        g[v + n_l].push_back(u);
      else
        g[u].push_back(v + n_l);
    }

    std::vector<int> queue;

    for (int i = 0; i < n_l; i++)
      if (vis[i]) queue.push_back(i);

    for (int qq = 0; qq < (int)queue.size(); ++qq) {
      int u = queue[qq];
      for (auto &v : g[u]) {
        if (vis[v]) continue;
        vis[v] = 1;
        queue.push_back(v);
      }
    }

    for (auto &pack : M) {
      auto &u = pack[0];
      auto &v = pack[0];
      if (vis[v])
        cover.push_back(v);
      else
        cover.push_back(u);
    }
    return cover;
  }
};
