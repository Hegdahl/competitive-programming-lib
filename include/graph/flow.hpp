#pragma once

#include <algorithm>
#include <array>
#include <limits>
#include <local_assertion.hpp>
#include <numeric>
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
  // hi! I have :floshed: idea for posets we do
  // Poset poset(n); poset(i) < poset(j); makes edge i -> j,
  // and all other implied edges

  // orz idea :happiness:
  // ur writing too wide lol i can't scroll that far for some reason
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

template <class FlowType = long long,
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
  std::vector<std::vector<int>> graph, bfs_graph, bfs_graph_rev;
  std::vector<int> vis;

  template <class... Args>
  NodeReference operator()(Args &&...args) {
    return NodeReference{
        *this, node_name_translator.to_id(std::forward<Args>(args)...)};
  }

  template <class... Args>
  std::vector<std::tuple<int, int, FlowType, FlowType>> get_edge_group(int x, int y) {
    auto [from, is_to] = node_name_translator.get_edge_group_descriptor(x, y);

    std::vector<std::tuple<int, int, FlowType, FlowType>> res;
    for (int i : from) {
      for (int eid : graph[i]) {
        Edge &edge = edge_list[eid];
        int j = edge.to;
        if (!is_to(j)) continue;

        FlowType current_flow = edge_list[eid^1].cap;
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

    vis.assign(graph.size(), 0);

    int steps = 0;

    FlowType ans = 0;
    FlowType to_add;
    while ((to_add = blocking_flow(source_id, sink_id))) ans += to_add, ++steps;

    return ans;
  }

 private:
  void expand_graph(int new_size) {
    while ((int)graph.size() < new_size) graph.emplace_back();
  }

  /**
   * @brief Compute flow in a BFS graph
   * until there are no unsaturated
   * paths left inside the bfs graph.
   *
   * @param source_id The id of the node to start from
   * @param sink_id   The id of the node to end at
   * @return The maximum flow without using edges
   *         outside the BFS graph.
   */
  FlowType blocking_flow(int source_id, int sink_id) {
    init_bfs_graph(source_id, sink_id);

    int steps = 0;

    FlowType ans = 0;
    FlowType to_add;
    while ((to_add =
                dfs(source_id, sink_id, std::numeric_limits<FlowType>::max())))
      ans += to_add, ++steps;

    return ans;
  }

  void init_bfs_graph(int source_id, int sink_id) {
    static std::vector<int> queue;

    if (bfs_graph.size() < graph.size()) {
      bfs_graph.assign(graph.size(), {});
      bfs_graph_rev.assign(graph.size(), {});
    }

    for (int i : queue) {
      bfs_graph[i].clear();
      bfs_graph_rev[i].clear();
    }

    queue.assign({source_id});

    vis.assign(graph.size(), 0);
    vis[source_id] = 1;

    for (int qq = 0; qq < (int)queue.size(); ++qq) {
      int i = queue[qq];
      if (i == sink_id)
        break;

      for (int eid : graph[i]) {
        Edge &edge = edge_list[eid];
        if (!edge.cap) continue;

        int j = edge.to;

        if (!vis[j]) {
          vis[j] = vis[i] + 1;
          queue.push_back(j);
        }
        if (vis[i] + 1 != vis[j]) continue;

        bfs_graph_rev[j].push_back(eid);
        //bfs_graph[i].push_back(eid);
      }
    }

    for (int i : queue)
      vis[i] = 0;
    vis[sink_id] = 1;

    for (int qq = (int)queue.size()-1; qq >= 0; --qq) {
      int j = queue[qq];

      if (!vis[j])
        continue;

      for (int eid : bfs_graph_rev[j]) {
        int i = edge_list[eid^1].to;
        vis[i] = 1;
        bfs_graph[i].push_back(eid);
      }
    }

    for (int i : queue)
      vis[i] = 0;
  }

  FlowType dfs(int current, int target, FlowType prefix_cap) {
    if (current == target) return prefix_cap;

    while (bfs_graph[current].size()) {
      Edge &edge = edge_list[bfs_graph[current].back()];
      int next = edge.to;

      if (vis[next]) {
        bfs_graph[current].pop_back();
        continue;
      }

      if (!edge.cap) {
        bfs_graph[current].pop_back();
        continue;
      }

      FlowType here = dfs(next, target, std::min(prefix_cap, edge.cap));
      if (!here) {
        bfs_graph[current].pop_back();
        continue;
      }

      edge.cap -= here;
      edge_list[bfs_graph[current].back()^1].cap += here;

      if (!edge.cap) {
        bfs_graph[current].pop_back();
        if (bfs_graph[current].size() == 0)
          vis[current] = 1;
      }

      return here;
    }

    vis[current] = 1;
    return 0;
  }
};