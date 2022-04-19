#pragma once

#include <local_assertion.hpp>

#include <array>
#include <vector>

/**
 * @brief Translate (layer_number, index_in_layer)
 * to node index for fixed size layers.
 */
struct LayeredNodeNameTranslator {
  const std::vector<int> layer_sizes;
  const std::vector<int> layer_sizes_prefix_sum;

  template<class... Args>
  LayeredNodeNameTranslator(Args &...args)
      : layer_sizes(std::forward<Args>(args)...) {
    layer_sizes_prefix_sum.resize(layer_sizes.size()+1);
    for (int i = 0; i < (int)layer_sizes.size(); ++i)
      layer_sizes_prefix_sum[i+1] =
        layer_sizes_prefix_sum[i] + layer_sizes[i];
  }

  /**
   * @brief Find node name
   * 
   * @param layer_number   which layer to find a node in
   * @param index_in_layer 0 for first,
   *                       layer_size[layer_number]-1 for last
   * @return node index
   */
  int operator()(int layer_number, int index_in_layer) {
    LOCAL_ASSERT_IN_RANGE(layer_number, 0, (int)layer_sizes.size()-1);
    LOCAL_ASSERT_IN_RANGE(index_in_layer, 0, layer_sizes[layer_number]-1);
    return layer_sizes_prefix_sum[layer_number] + index_in_layer;
  }

  /**
   * @brief Assume layer zero if no layer is given
   * 
   * @param index_in_layer index in layer zero
   * @return node index
   */
  int operator()(int index_in_layer) {
    return (*this)(0, index_in_layer);
  }
};

template<class FlowType = long long, class NodeNameTranslator = LayeredNodeNameTranslator>
class Flow {

 public:

  NodeNameTranslator node_name_translator;
  using NodeNameTranslator::NodeNameTranslator;

};
