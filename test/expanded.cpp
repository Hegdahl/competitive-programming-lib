/**
 * problem: F. Making It Bipartite
 * url:     https://codeforces.com/contest/1630/problem/F
 */
// <graph/flow.hpp>
// <local_assertion.hpp>
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
               << '[' << low << ", " << high << ']')// </local_assertion.hpp>
// <utils.hpp>
#include <random>
#include <chrono>
std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());
struct r_hash {
  static uint64_t splitmix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
    x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
    return x ^ (x >> 31);
  }
  template<class T>
  uint64_t operator()(const T &x) const;
  uint64_t operator()(uint64_t x) const {
    static const uint64_t FIXED_RANDOM = rng();
    return splitmix64(x + FIXED_RANDOM);
  }
};
template<class...Ts>
struct overloaded : Ts... { using Ts::operator()...; };
template<class... Ts>
overloaded(Ts...) -> overloaded<Ts...>;
// </utils.hpp>
// <container/any_map.hpp>
#include <map>
#include <tuple>
#include <utility>
template<class T>
class AnyMap {
  static int instance_count;
  const int instance_id_;
  T default_value_;
 public:
  AnyMap(T default_value = T())
      : instance_id_(instance_count++), default_value_(std::move(default_value)) {}
  template <class... Args>
  auto &operator()(Args &&...args) {
    using Tuple =
        decltype(std::make_tuple(instance_id_, std::forward<Args>(args)...));
    static std::map<Tuple, T> map;
    auto [it, placed] = map.emplace(
        Tuple{instance_id_, std::forward<Args>(args)...}, default_value_);
    return it->second;
  }
};
template<class T>
int AnyMap<T>::instance_count = 0;// </container/any_map.hpp>
#include <algorithm>
#include <array>
#include <limits>
#include <numeric>
#include <string>
#include <tuple>
#include <vector>
/**
 * NodeNameTranslator interface
 *
 * needed:
 * to_id(...):
 *  should translate some kind of name to a node id.
 *
 * optional:
 * from_id(int):
 *  should translate node id to name.
 * get_edge_group_descriptor(...):
 *  should return a list of node ids representing
 *  starting endpoints in the edge group,
 *  and a function object determining
 *  returing true given a valid node id for
 *  an edge endpoint in the edge group.
 */
/**
 * @brief Translate (layer_number, index_in_layer)
 * to node index for fixed size layers.
 */
struct LayeredNodeNameTranslator {
  const std::vector<int> layer_sizes;
  const std::vector<int> layer_sizes_prefix_sum;
  template <class... Args>
  LayeredNodeNameTranslator(Args &&...args)
      : layer_sizes{std::forward<Args>(args)...},
        layer_sizes_prefix_sum(get_prefix_sum(layer_sizes)) {}
  static std::vector<int> get_prefix_sum(const std::vector<int> &a) {
    std::vector<int> res(a.size() + 1);
    for (int i = 0; i < (int)a.size(); ++i) res[i + 1] = res[i] + a[i];
    return res;
  }
  /**
   * @brief Find node name
   *
   * @param layer_number   which layer to find a node in
   * @param index_in_layer 0 for first,
   *                       layer_size[layer_number]-1 for last
   * @return node index
   */
  int to_id(int layer_number, int index_in_layer) {
    LOCAL_ASSERT_IN_RANGE(layer_number, 0, (int)layer_sizes.size() - 1);
    LOCAL_ASSERT_IN_RANGE(index_in_layer, 0, layer_sizes[layer_number] - 1);
    return layer_sizes_prefix_sum[layer_number] + index_in_layer;
  }
  /**
   * @brief Assume layer zero if no layer is given
   *
   * @param index_in_layer index in layer zero
   * @return node index
   */
  int to_id(int index_in_layer) { return to_id(0, index_in_layer); }
  std::array<int, 2> from_id(int id) {
    LOCAL_ASSERT_IN_RANGE(id, 0, layer_sizes_prefix_sum.back() - 1);
    auto it = std::upper_bound(layer_sizes_prefix_sum.begin(),
                               layer_sizes_prefix_sum.end(), id);
    LOCAL_ASSERT(it != layer_sizes_prefix_sum.begin(),
                 "oops, something is wrong with the internal binary search.");
    int layer_number = int(it - layer_sizes_prefix_sum.begin()) - 1;
    int index_in_layer = id - layer_sizes_prefix_sum[layer_number];
    return {layer_number, index_in_layer};
  }
  auto get_edge_group_descriptor(int from_layer, int to_layer) {
    std::vector<int> from_list(layer_sizes[from_layer]);
    std::iota(from_list.begin(), from_list.end(),
              layer_sizes_prefix_sum[from_layer]);
    auto is_in_to_layer =
        [begin = layer_sizes_prefix_sum[to_layer],
         end = layer_sizes_prefix_sum[to_layer + 1]](int index) {
          return begin <= index && index < end;
        };
    return std::make_pair(from_list, is_in_to_layer);
  }
};
struct AnyNodeNameTranslator {
  AnyMap<int> mp;
  int next_node_id;
  AnyNodeNameTranslator() : mp(-1), next_node_id(0) {}
  template <class... Args>
  int to_id(Args &&...args) {
    int &v = mp(std::forward<Args>(args)...);
    if (v == -1) v = next_node_id++;
    return v;
  }
};
template <class FlowType = long long, bool shuffle = false,
          class NodeNameTranslator = LayeredNodeNameTranslator>
class Flow {
  struct EdgeReference {
    Flow &flow;
    int i, j;
    void operator=(FlowType w) { flow.make_edge(i, j, w); }
    FlowType max_flow() { return flow.max_flow(i, j); }
  };
  struct NodeReference {
    Flow &flow;
    int i;
    template <class... Args>
    EdgeReference operator()(Args &&...args) {
      return {flow, i,
              flow.node_name_translator.to_id(std::forward<Args>(args)...)};
    }
  };
 public:
  NodeNameTranslator node_name_translator;
  template <class... Args>
  Flow(Args &&...args) : node_name_translator(std::forward<Args>(args)...) {}
  struct Edge {
    int to;
    FlowType cap;
  };
  std::vector<Edge> edge_list;
  std::vector<std::basic_string<int>> graph;
  std::vector<int> bfs_graph, bfs_graph_outdegree, bfs_graph_start, queue,
      used_edge;
  std::vector<int> vis, dist;
  template <class... Args>
  NodeReference operator()(Args &&...args) {
    return NodeReference{
        *this, node_name_translator.to_id(std::forward<Args>(args)...)};
  }
  template <class... Args>
  std::vector<std::tuple<int, int, FlowType, FlowType>> get_edge_group(int x,
                                                                       int y) {
    auto [from, is_to] = node_name_translator.get_edge_group_descriptor(x, y);
    std::vector<std::tuple<int, int, FlowType, FlowType>> res;
    for (int i : from) {
      for (int eid : graph[i]) {
        Edge &edge = edge_list[eid];
        int j = edge.to;
        if (!is_to(j)) continue;
        FlowType current_flow = edge_list[eid ^ 1].cap;
        FlowType capacity = edge.cap + current_flow;
        res.emplace_back(i, j, current_flow, capacity);
      }
    }
    return res;
  }
  /**
   * @brief Creates an edge from `source_id` to `sink_id`
   * with capacity equal to `weight`.
   *
   * To use node names instead, use this:
   * `flow_instance(source node name)(sink node name) = weight`
   * which uses this method internally.
   *
   * @param from_id The starting node of the edge
   * @param to_id   The ending node of the edge
   * @param weight  The capacity of the edge
   */
  void make_edge(int from_id, int to_id, FlowType weight) {
    LOCAL_ASSERT(from_id != to_id, "Self-loops are not supported.")
    expand_graph(std::max(from_id, to_id) + 1);
    graph[from_id].push_back(edge_list.size());
    edge_list.emplace_back(Edge{to_id, weight});
    graph[to_id].push_back(edge_list.size());
    edge_list.emplace_back(Edge{from_id, 0});
  }
  /**
   * @brief Compute the maximum flow from
   * `source_id` to `sink_id`.
   *
   * To use node names instead, use this:
   * `flow_instance(source node name)(sink node name).max_flow()`
   * which uses this method internally.
   *
   * @param source_id The id of the node to start from
   * @param sink_id   The id of the node to end at
   * @return The maximum flow
   */
  FlowType max_flow(int source_id, int sink_id) {
    expand_graph(std::max(source_id, sink_id) + 1);
    source_id_ = source_id;
    sink_id_ = sink_id;
    if constexpr (shuffle)
      for (auto &adj_list : graph)
        std::shuffle(adj_list.begin(), adj_list.end(), rng);
    FlowType ans = 0;
    FlowType to_add;
    while ((to_add = blocking_flow())) ans += to_add;
    return ans;
  }
 private:
  int source_id_, sink_id_;
  void expand_graph(int new_size) {
    while ((int)graph.size() < new_size) {
      graph.emplace_back();
      bfs_graph_outdegree.emplace_back();
      bfs_graph_start.emplace_back();
      vis.emplace_back();
      dist.emplace_back();
    }
  }
  /**
   * @brief Compute flow in a BFS graph
   * until there are no unsaturated
   * paths left inside the bfs graph.
   *
   * @return The maximum flow without using edges
   *         outside the BFS graph.
   */
  FlowType blocking_flow() {
    init_bfs_graph();
    FlowType ans = 0;
    FlowType to_add;
    while ((to_add =
                dfs(source_id_, std::numeric_limits<FlowType>::max())))
      ans += to_add;
    return ans;
  }
  void init_bfs_graph() {
    reset_bfs_graph();
    forward_bfs();
    backward_bfs();
    load_edges();
  }
  void reset_bfs_graph() {
    for (int i : queue) {
      vis[i] = 0;
      dist[i] = 0;
      bfs_graph_outdegree[i] = 0;
    }
    queue.clear();
  }
  void forward_bfs() {
    queue.push_back(source_id_);
    vis[source_id_] = 1;
    for (int qq = 0; qq < (int)queue.size(); ++qq) {
      int i = queue[qq];
      if (i == sink_id_) break;
      for (int eid : graph[i]) {
        Edge &edge = edge_list[eid];
        if (!edge.cap) continue;
        int j = edge.to;
        if (vis[j]) continue;
        vis[j] = 1;
        dist[j] = dist[i] + 1;
        queue.push_back(j);
      }
    }
  }
  void backward_bfs() {
    vis[sink_id_] = 0;
    bfs_graph_start[queue.back()] = 0;
    for (int qq = (int)queue.size() - 1; qq >= 0; --qq) {
      int j = queue[qq];
      if (qq + 1 < (int)queue.size()) {
        int nxt = queue[qq + 1];
        bfs_graph_start[j] = bfs_graph_start[nxt] + bfs_graph_outdegree[nxt];
        bfs_graph_outdegree[nxt] = 0;
      }
      if (vis[j]) continue;
      for (int eid : graph[j]) {
        if (!edge_list[eid ^ 1].cap) continue;
        int i = edge_list[eid].to;
        if (dist[i] + 1 != dist[j]) continue;
        vis[i] = 0;
        ++bfs_graph_outdegree[i];
        used_edge.push_back(eid);
      }
    }
  }
  void load_edges() {
    bfs_graph.resize(bfs_graph_start[queue.front()] +
                     bfs_graph_outdegree[queue.front()]);
    bfs_graph_outdegree[queue.front()] = 0;
    while (used_edge.size()) {
      int rev_eid = used_edge.back();
      used_edge.pop_back();
      int i = edge_list[rev_eid].to;
      bfs_graph[bfs_graph_start[i] + bfs_graph_outdegree[i]++] = rev_eid ^ 1;
    }
  }
  FlowType dfs(int current, FlowType prefix_cap) {
    if (current == sink_id_) return prefix_cap;
    while (bfs_graph_outdegree[current]) {
      int eid = bfs_graph[bfs_graph_start[current] +
                          bfs_graph_outdegree[current] - 1];
      Edge &edge = edge_list[eid];
      int next = edge.to;
      if (vis[next]) {
        --bfs_graph_outdegree[current];
        continue;
      }
      if (!edge.cap) {
        --bfs_graph_outdegree[current];
        continue;
      }
      FlowType here = dfs(next, std::min(prefix_cap, edge.cap));
      if (!here) {
        --bfs_graph_outdegree[current];
        continue;
      }
      edge.cap -= here;
      edge_list[eid ^ 1].cap += here;
      if (!edge.cap) {
        if (!--bfs_graph_outdegree[current]) vis[current] = 1;
      }
      return here;
    }
    vis[current] = 1;
    return 0;
  }
};// </graph/flow.hpp>
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
