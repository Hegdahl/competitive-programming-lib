#pragma once

#include <map>
#include <tuple>
#include <utility>

template <class T>
class AnyMap {
  static int instance_count;
  const int instance_id_;
  T default_value_;

 public:
  AnyMap(T default_value = T())
      : instance_id_(instance_count++),
        default_value_(std::move(default_value)) {}

  template <class... Args>
  auto &operator()(Args &&...args) {
    using Tuple =
        decltype(std::make_tuple(instance_id_, std::forward<Args>(args)...));
    static std::map<Tuple, T> map;

    auto it = map.emplace(Tuple{instance_id_, std::forward<Args>(args)...},
                          default_value_)
                  .first;
    return it->second;
  }
};

template <class T>
int AnyMap<T>::instance_count = 0;
