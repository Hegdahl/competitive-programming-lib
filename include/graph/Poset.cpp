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
    Flow<int> internal_flow;
    int n;
    Bipar
    Poset(int n) : Bipartite_matcher(1, n, n, 1) { init(); }
    void init(){
        //for(int i=0;i<n_l;i++) internal_flow(0, 0)(1, i) = 1;
        //for(int i=0;i<n_r;i++) internal_flow(2, i)(3, 0) = 1;
    }

    // how to tell the difference between
    // flow(i_want_this_many_layers, i_want_each_layer_to_be_this_big)
    // and
    // flow(layer_size_0, layer_size_1)
    // flow(layer_size_0, layer_size_1, layer_size_2)
    // i tried just forwarding to vector constructor
    // so that the second type uses {}, but a consequence
    // of the forwarding is turing {} into () :angeryboye:
    // maybe the first type isn't useful and we can just always use 2nd type
    //good question :thinkies:
    void add_edge(int u, int v){
        //internal_flow(1, u)(2, v) = 1;
    }
    int compute_matching_size(){
        //return internal_flow(0, 0)(3, 0).max_flow();
    }
    std::vector<std::array<int, 2>> get_matching(){
        int matching_size = compute_matching_size();
        std::vector<std::array<int, 2>> matchings(matching_size);
        for(const auto &[u, v, flow, cap] : internal_flow.get_edge_group(1, 2)){
            if(flow > 0) matching.push_back(array<int, 2>({u, v})); //updating names later
        }
        return matching;
    }
    vector<int> min_vertex_cover(){
        auto M = get_matching();
        std::vector<int> cover(M.size());
        std::vector<uint8_t> vis(n_l + n_r);
        for(auto &[u, v] : M) vis[u] = 1;
        std::vector<std::vector<int>> g(n_l + n_r);
        for(auto &[u, v, flow, cap] : internal_flow.get_edge_group(1, 2)){
            if(flow) g[v + n_l].push_back(u);
            else g[u].push_back(v + n_l);
        }
        std::queue<int> bfs;
        for(int i=0;i<n_l;i++) if(vis[i]) bfs.push(i);
        for(int i=0;i<n_l;i++){
            if(is_source[i]) bfs.push(i), vis[i] = 1;
        }
        while(bfs.size()){
            int u = bfs.front();
            bfs.pop();
            for(auto &v : g[u]){
                if(vis[v]) continue;
                vis[v] = 1;
                bfs.push(v);
            }
        }
        for(auto &[u, v] : M){
            if(vis[v] == 1){
                cover.push_back(v);
            }
            else{
                cover.push_back(u);
            }
        }
        return cover;
    }

    // there is actually 1 feature of java i'm slightly missing now:
    // being able to initialize const members in the body of the constructor :da:
    
    //I see smart
    //so that reduces mem from 8 -> 4 yea
    // i have funny idea
    // what if we don't store edge.to in the edges, but instead edge.toggle
    // then for some i,j  the edge goes from i to i^edge.toggle
    //can you store edge refernce instead?
    // can't have references in a vector, but we can have indices (edge ids) to edges
    // I can't do structured binding iteration then
    // yes you can, the get_edge_group function doesn't return real edges.
    // just use it as you already do
    // :O didn't think that was possible
};

// should we have custom assertion system that is only enabled locally?
// 
//seems pretty easy I think, why not
// basically it's just this:
// you can assume any node exists
// then just (x, y, z, w) will refer to some node
// with any number of indices
//wait why would you do that?
//I think it makes more sense {layer number, index_in_layer}
//tho your idea is superset of my idea
// ^ but worse constant factor. but i think the constant won't matter
// since flow is slow compared to lookups in map of arrays
//log * num_layer is still bad imo 
//maybe y

// it's log * number of dimensions, not log * number of layers
//yea sorry thats what I meant, thats still pretty bad imo
//you can do in linear in number of nodes, if you assume it to behave like nested vectors
//i.e. nested layers, for 3d (x, y, z) find all nodes mapped to by x which, in some prefix [1, k]
//then it reduces to k 2d cases
//but this won't allow you to use x = 0, 1, 1e9
//its essentially indexing by radix sort on (x, y, z)

// i think i can abstract out this translation thing
// to have the best of both worlds
//ok I'm watching :blob_popcorn:

// just assume it works like you thought it would initially :thumbsup:
//ok trust
// with this mapping idea we actually don't need to predefine
// size of graph

//take size of layer = maximum index + 1(if 0 indexed)?
//do you use 1 based indexing?
// ^ up to user
// if i open and close files, do they open and close for you too?
//no, it doesn't
//ah I see what you mean, we can have them built dynamically
// ^^
//I will need to see your impl a bit, to see how to do, that, I will edit parts accordingly
