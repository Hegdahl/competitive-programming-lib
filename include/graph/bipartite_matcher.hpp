#include <local_assertion.hpp>
#include <graph/flow.hpp>
#include<array>
#include<vector>
/**
 * @brief computes the following :
 *    maximum matching on bipartite graph
 *    minimum vertex cover of bipartite graph
*/
struct Bipartite_matcher{
    Flow<int> internal_flow;
    int n_l, n_r;
    Bipartite_matcher(int n) : n_l(n), n_r(n), internal_flow(1, n, n, 1) {}
    Bipartite_matcher(int n_l, int n_r) : n_l(n_l), n_r(n_r), internal_flow(1, n_l, n_r, 1) {}
    void add_edge(int u, int v){
        internal_flow(1, u)(2, v) = 1; 
    }
    int compute_matching_size(){
        //TO DO add edges from source and sink
        return internal_flow.max_flow();
    }
    std::vector<std::array<int, 2>> get_matching(){
        std::vector<std::array<int, 2>> matchings;
    }
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
