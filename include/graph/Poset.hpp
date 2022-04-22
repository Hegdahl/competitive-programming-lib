#pragma once

#include <graph/bipartite_matcher.hpp>
#include <graph/flow.hpp>
#include <local_assertion.hpp>

#include <array>
#include <vector>

/**
 * @brief Poset container
 * computes min number/max size chains/antichains
*/
struct Poset{
    int n;
    struct Poset_node_reference{
        Poset &poset;
        int u;
        Poset_node_reference(Poset &poset_, int u) : poset(poset_), u(u) {}
        void operator<(Poset_node_reference &o){
            poset.add_edge(u, o.u);
        }
    };
    Poset_node_reference operator()(int x){
        return Poset_node_reference(*this, x);
    }
    void add_edge(int u,int v){

    }
    Poset(int n) : BipartiteMatcher(n, n) { init(); }
};
