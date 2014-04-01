/**
 * @file forest_fire_model.cpp
 * Implement a simple version of the forest fire model.
 * @author David F. Gleich
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include <vector>
#include <algorithm>
#include <iterator>
#include <set>

#include "argstream.h"

#include "sparfun_util.h"


typedef size_t vertex_type;

/** 
 * A class to control a forest fire model of a graph.
 *
 * Usage:
 * 
 *  forest_fire_model M(0.1); // burning probability
 *  // initialize the seed graph with a clique of two vertices.
 *  M.init(2, 10); // also reserve space for 10 extra vertices.
 *  // now add a series of vertices
 *  for (int i=0; i<10; ++i) { M.new_vertex(); }
 *  // get the edges
 *  for (size_t i=0; i<M.nverts; ++i) {
 *    size_t neigh = M.degree(i);
 *    while(neigh-->0) { size_t j = J.neighbor(i,neigh);  // (i,j) is an edge }
 *  }
 */
class forest_fire_model
{
    
public:
    double alpha;
    vertex_type nverts; 
    vertex_type finished_vert;
    size_t nedges;
    std::vector< std::vector<vertex_type> > graph;
    std::vector< vertex_type > parent;
    forest_fire_model(double alpha_)
    : alpha(alpha_), nverts(0), finished_vert(0)
    {}
    
    /** Reset the data structure */
    void reset(size_t newsize) { 
        for (size_t i = 0; i<nverts; ++i) {
            graph[i].resize(0);
        }
        graph.reserve(newsize);
        graph.resize(0);
        parent.reserve(newsize);
        parent.resize(0);
        finished_vert = 0;
        nedges = 0;
    }
        
    vertex_type alloc_vertex(vertex_type p) {
        vertex_type rval = nverts;
        nverts+=1;
        assert(p < nverts); // check here, because p could be it's own parent
        parent.push_back(p);
        assert(parent.size() == nverts);
        // use push_back if we exceed our capacity
        if (nverts > graph.capacity()) {
            graph.push_back(std::vector<vertex_type>());
        } else {
            graph.resize(nverts);
        }
        assert(graph.size() == nverts);
        return rval;
    }
    
    void add_edge(vertex_type u, vertex_type v) {
        assert(u < v);
        if (v > finished_vert) {
            // then we do not have to search for edges
            graph[v].push_back(u);
            graph[u].push_back(v);
            nedges ++;
        } else {
            // here we have to search for edges
            if (!is_edge(u,v)) {
                graph[v].push_back(u);
                graph[u].push_back(v);
                nedges ++;
            }
        }
    }
    
    /** Initialize the model with a clique graph 
     * @param num_new_verts the number of vertices added 
     */
    void init(size_t clique_size, size_t num_new_verts) {
        reset(clique_size + num_new_verts);
        for (size_t v = 0; v < clique_size; ++v) {
            vertex_type i = alloc_vertex(v);
            for (size_t u=0; u < v; ++u) {
                add_edge(u,v);
            }
            finished_vert = i;
        }
    }
    
    size_t degree(vertex_type i) {
        return graph[i].size();
    }
    
    /** Return the jth neighbor of vertex i */
    vertex_type neighbor(vertex_type i, size_t j) {
        assert(j < degree(i));
        return graph[i][j];
    }
    
    bool is_edge(vertex_type i, vertex_type j) {
        if (degree(i) < degree(j)) {
            // search in i
            for (size_t k = 0; k<degree(i); k++) {
                if (neighbor(i,k)==j) {
                    return true;
                }
            }
        } else {
            // search in j
            for (size_t k = 0; k<degree(j); k++) {
                if (neighbor(j,k)==i) {
                    return true;
                }
            }
        }
        return false;
    }
    
    void burn(vertex_type start, std::vector<vertex_type>& list) {
        list.clear();
        typedef std::set<vertex_type> set_type;
        set_type burning; // the current phase of burning
        set_type burned; // includes all of burning
        set_type toburn; // the next phase of burning
        std::vector<vertex_type> alive;
        
        toburn.insert(start);
        burned.insert(start);
        list.push_back(start);
        
        while(toburn.size() > 0) {
            burning = toburn;
            toburn.clear();
            for (set_type::const_iterator burnit=burning.begin(); 
                    burnit != burning.end(); ++burnit) {
                vertex_type v = *burnit;
                //assert(burned.count(v) == 0);
                //burned.insert(v);
                
                // insert the unburned neighbors into a list
                alive.clear();
                for (size_t k = 0; k < degree(v); ++k) {
                    vertex_type u = neighbor(v,k);
                    if (burned.count(u) == 0) {
                        alive.push_back(u);
                    }
                }
                
                // randomly shufle the list
                std::random_shuffle( alive.begin(), alive.end() );
                size_t nlinks = std::min(sf_randgeo(1.0-alpha)-1, alive.size());
                for (size_t k = 0; k < nlinks; ++k) {
                    toburn.insert(alive[k]);
                    burned.insert(alive[k]);
                    list.push_back(alive[k]);
                }
            }
        }
    }
        
    
    /** Implement the routine for the copying model to add a vertex. */
    void add_vertex() {
        // pick a prototype.
        
        vertex_type proto = (vertex_type)sf_randint(0,nverts-1);
        
        std::vector<vertex_type> burned;
        burn(proto,burned);
        
        vertex_type i = alloc_vertex(proto);
        for (size_t k = 0; k < burned.size(); ++k) {
            add_edge(burned[k],i);
        }
        finished_vert = i; // we are done up to here!
    }
};

bool generate_graph(size_t initial_verts, double alpha, size_t total_verts,
    const char* outputfilename)
{
    forest_fire_model M(alpha);
    M.init(initial_verts, total_verts - initial_verts);
    for (size_t i=initial_verts; i<total_verts; ++i) { M.add_vertex(); }
    
    // count the number of edges
    size_t nedges = 0;
    for (vertex_type i=0; i<M.nverts; ++i) { nedges += M.degree(i); }
    
    // output the graph
    FILE* of = fopen(outputfilename,"wb");
    if (!of) { 
        fprintf(stderr, "error, cannot open file %s\n", outputfilename);
        return false;
    }
    fprintf(of, "%zi %zi %zi\n", total_verts, total_verts, nedges);
    for (vertex_type i=0; i<M.nverts; ++i) {
        std::vector< vertex_type > neighs(M.degree(i));
        for (size_t k=0; k<M.degree(i); ++k) {
            neighs[k] = M.neighbor(i,k);
        }
        std::sort(neighs.begin(), neighs.end());
        for (size_t k=0; k<M.degree(i); ++k) {
            fprintf(of, "%zi %zi 1\n", i, neighs[k]);
        }
    }
    fclose(of);
    return true;
}    

int main(int argc, char **argv) 
{
    // arguments: m (initial) , p, n (total)
    // output name
    
    size_t initial, total;
    double p;
    int seed=-1;
    bool force=false;
    std::string outputfilename;
    
    argstream::argstream as(argc, argv);
    as >> argstream::parameter('m',"initial",initial,
        "The number of initial vertices.", true);
        
    as >> argstream::parameter('n',"total",total,
        "The number of total vertices.", true);
        
    as >> argstream::parameter('p',"prob",p,
        "The burning probability.", true);
    
    as >> argstream::parameter('s',"seed",seed,
        "The random seed (-1 means seed from time) ", false);
    
    as >> argstream::parameter('o',"output",outputfilename,
        "The output filename", true);
        
    as >> argstream::option('f',"force",force,
        "Overwrite the file even if it exists.");
        
    as >> argstream::help();
    as.defaultErrorHandling();
    
    if (!force) {
        if (sf_file_exists(outputfilename.c_str())) {
            fprintf(stderr,
                "error: output file %s already exists (-f to overwrite)\n",
                outputfilename.c_str());
            return (-1);
        }
    }
    
    if (seed == -1) {
        sf_timeseed();
    } else {
        sf_srand((unsigned int)seed);
    }
    
    if (total < initial) {
        total = initial;
    }
    
    if (p < 0.0 || p > 1.0) {
        fprintf(stderr, "error p=%f is not in [0,1]\n", p);
        return (-1);
    }
    
    if (!generate_graph(initial, p, total, outputfilename.c_str())) {
        return (1);
    }
}



/*
    
    forest_fire_model M(alpha);
    M.init(init_verts,num_verts);
    for (size_t i=0; i<num_verts; ++i) {
        M.add_vertex();
    }
    
    // count the number of edges
    size_t en = 0;
    for (vertex_type i=0; i<M.nverts; ++i) { en += M.degree(i); }
    
    // copy data out
    dims[0] = en;
    dims[1] = 2;
    plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
    double *edata_src = (double*)mxGetData(plhs[0]);
    double *edata_dst = edata_src + en;
    
    en = 0;
    for (vertex_type i=0; i<M.nverts; ++i) {
        for (size_t k=0; k<M.degree(i); ++k) {
            edata_src[en] = (double)(i+1);
            edata_dst[en] = (double)(M.neighbor(i,k)+1);
            en++;
        }
        if (parray) {
            parray[i] = (double)(M.parent[i]+1);
        }
    }
}
*/
