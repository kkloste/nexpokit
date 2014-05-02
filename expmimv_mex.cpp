/**
 * @file expmimv_mex.cpp
 * @author Kyle Kloster, David F. Gleich
 */

/**
 * This file implements a gauss-southwell type method for the truncated
 * taylor series approximation for a column of the matrix exponential
 */



#ifndef __APPLE__
#define __STDC_UTF_16__ 1
#endif 

#include <sparsehash/dense_hash_map>


#include <vector>
#include <queue>
#include <utility> // for pair sorting
#include <assert.h>
#include <limits>
#include <algorithm>
#include <math.h>

#include "taydeg.hpp"

#include "mex.h"

#define DEBUGPRINT(x) do { if (debugflag) { \
mexPrintf x; mexEvalString("drawnow"); } \
} while (0)

int debugflag = 0;
    typedef google::dense_hash_map<mwIndex,double> map_type;
    
struct sparsevec {
    typedef google::dense_hash_map<mwIndex,double> map_type;
    map_type map;
    sparsevec()  {
        map.set_empty_key((mwIndex)(-1));
    }
    /** Get an element and provide a default value when it doesn't exist
     * This command does not insert the element into the vector
     */
    double get(mwIndex index, double default_value=0.0) {
        map_type::iterator it = map.find(index);
        if (it == map.end()) {
            return default_value;
        } else {
            return it->second;
        }
    }
    
    /** Compute the sum of all the elements
     * Implements compensated summation
     */
    double sum() {
        double s=0.;
        for (map_type::iterator it=map.begin(),itend=map.end();it!=itend;++it) {
            s += it->second;
        }
        return s;
    }
    
    /** Compute the max of the element values
     * This operation returns the first element if the vector is empty.
     */
    mwIndex max_index() {
        mwIndex index=0;
        double maxval=std::numeric_limits<double>::min();
        for (map_type::iterator it=map.begin(),itend=map.end();it!=itend;++it) {
            if (it->second>maxval) { maxval = it->second; index = it->first; }
        }
        return index;
    }
};

mwIndex heap_up(mwIndex j, mwSize n, mwIndex* T, double* d) {
    while (1) {
        if (j==0) { break; } /* the element is at the top */
        mwIndex j2 = (j-1)/2;
        double valj2 = d[j2];
        if (valj2 < d[j]) {
            break; /* the parent is larger, so stop */
        } else {
            /* the parent is larger, so swap */
            mwIndex heapj = T[j];
            mwIndex heapj2 = T[j2];
            d[j2] = d[j];
            d[j] = valj2;
            T[j2] = heapj;// L[heapj] = j2;
            T[j] = heapj2;// L[heapj2] = j;
            j = j2;
        }
    }
    return j;
}

/** Move an element down the heap until it hits the correct place.
 * i.e. it is bigger than its parent (small entries rise to the top)
 */
mwIndex heap_down(mwIndex k, mwSize n, mwIndex* T, double* d) {
    mwIndex heapk = T[k];
    double valk = d[k];
    while (1) {
        mwIndex i=2*(k+1)-1;
        if (i>=n) { break; } /* end of heap */
        if (i<n-1) {
            /* pick between children (unneeded if i==heap_size-1) */
            //            mwIndex left=T[i];
            //            mwIndex right=T[i+1];
            if (d[i+1] < d[i]) {
                i=i+1; /* pick the smaller child */
            }
        }
        if (d[k] < d[i]) {
            /* k is smaller than both children, so end */
            break;
        } else {
            T[k] = T[i];// L[T[i]]=k;
            T[i] = heapk;// L[heapk] = i;
            d[k] = d[i];
            d[i] = valk;
            k=i;
        }
    }
    return k;
}

struct sparserow {
    mwSize n, m;
    mwIndex *ai;
    mwIndex *aj;
    double *a;
};


/**
 * Returns the degree of node u in sparse graph s
 */
mwIndex sr_degree(sparserow *s, mwIndex u) {
    return (s->ai[u+1] - s->ai[u]);
}


/**
 * @param G - sparse matrix values (CSC)
 * @param set - the vector v to approximate exp(G)v (0 <= c < n)
 * @param y - the output vector (length n)
 * @param t - scaling parameter (0< t < inf)
 * @param eps - the stopping tolerance (0 < tol < Inf)
 * @param npushes - the number of operations (an output)
 */

void expm_svec(sparserow* G, std::vector<mwIndex>& set, sparsevec& y,
               const double t, const double tol,
               mwIndex maxnnz, double* npushes)
{
	double eps = tol;
    mwIndex N = taylordegree(t,eps);
    mwIndex n = G->n;
    
    DEBUGPRINT( ("Input n=%i N=%i tol=%f \n", n, N, eps) );
    
    // Initialize y = set
    for (size_t i=0; i<set.size(); ++i) {
        mwIndex ri = set[i];
        double rij = 1.;
        y.map[ri] = rij;
    }
    
    sparsevec dummy; // dummy temporarily holds A*y
    //    tree_map T;
    std::vector<mwIndex> Tvec(maxnnz,0);
    std::vector<double> dvec(maxnnz,0);
    mwIndex *T = &Tvec[0];
    double *d = &dvec[0];
    
    // HORNER'S RULE ITERATES
    for (mwIndex k=0; k <= N-1 ; k++){
        mwIndex hsize = 0;

        // (1) MAKE HEAP
        //      For each entry of y, check if we put it in the heap.
        for (map_type::iterator it=y.map.begin(),itend=y.map.end();it!=itend;++it) {
            mwIndex ind = it->first;
            double valind = it->second;
            
            // if the value is nonzero, we might add it to heap
            if (valind > 0) {
                // if the heap is full, we might have to delete
                if ( hsize >= maxnnz ){
                    double minval = d[0]; // y.map[T[0]];
                    if (valind > minval){
                        T[0] = ind; // replace root with new entry
                        d[0] = valind;
                        heap_down(0, hsize, T, d);
                        minval = d[0]; // y.map[T[0]];
                    }
                    else y.map[ind] = 0; // y[ind] so small it's excluded from heap
                }
                // if heap is not full, just place y[ind] at back of heap, then heap_up
                else {
                    T[hsize] = ind;
                    d[hsize] = valind;
                    hsize++;
                    heap_up(hsize-1, hsize, T, d);
                }
            }
        }// END FOR
        
        y.map.clear();
        
        // At this point, y[] = 0 entirely, and d[] contains
        // just the entries from y that survived the heap process,
        // i.e. the maxnnz largest entries that were in y.
        
        // (2) ACTUAL MATVEC
        for (mwIndex it = 0; it < hsize ; it++){
            mwIndex Tind = T[it];
            double Nk = (double)N-k;
            double valind = t*(d[it])/Nk;
            double update = valind/(double)sr_degree(G,Tind);
            *npushes += sr_degree(G,Tind);
            for (mwIndex nzi=G->ai[Tind]; nzi < G->ai[Tind+1]; ++nzi) {
                mwIndex v = G->aj[nzi];
                y.map[v] = y.get(v) + update;
            }
        }
        
        // add y += set
        for (size_t i=0; i<set.size(); ++i) {
            mwIndex ri = set[i];
            double rij = 1.;
            y.map[ri] = y.get(ri) + rij;
        }
        *npushes += set.size();
    }//terms of Taylor are complete
    return;
}

void copy_array_to_index_vector(const mxArray* v, std::vector<mwIndex>& vec)
{
    mxAssert(mxIsDouble(v), "array type is not double");
    size_t n = mxGetNumberOfElements(v);
    double *p = mxGetPr(v);
    
    vec.resize(n);
    
    for (size_t i=0; i<n; ++i) {
        double elem = p[i];
        mxAssert(elem >= 1, "Only positive integer elements allowed");
        vec[i] = (mwIndex)elem - 1;
    }
}


// USAGE
// [y npushes] = expmimv_mex(A,set,eps,t,maxnnz,debugflag)
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    if (nrhs < 2 || nrhs > 6) {
        mexErrMsgIdAndTxt("expmimv_mex:wrongNumberArguments",
                          "expmimv_mex needs two to six arguments, not %i", nrhs);
    }
    mxAssert(nlhs <= 2, "Too many output arguments");
    if (nrhs == 6) {
        debugflag = (int)mxGetScalar(prhs[5]);
    }
    DEBUGPRINT(("\n expmimv_mex: preprocessing start: \n"));
    
    const mxArray* mat = prhs[0];
    const mxArray* set = prhs[1];
    
    if ( mxIsSparse(mat) == false ){
        mexErrMsgIdAndTxt("expmimv_mex:wrongInputMatrix",
                          "expmimv_mex needs sparse input matrix");
    }
    if ( mxGetM(mat) != mxGetN(mat) ){
        mexErrMsgIdAndTxt("expmimv_mex:wrongInputMatrixDimensions",
                          "expmimv_mex needs square input matrix");
    }

    sparserow G;
    G.m = mxGetM(mat);
    G.n = mxGetN(mat);
    G.ai = mxGetJc(mat);
    G.aj = mxGetIr(mat);
    G.a = mxGetPr(mat);

    double* npushes = 0;
    if (nlhs > 1){
        plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
        npushes = mxGetPr(plhs[1]);
    }
    
    double eps = pow(10,-5);
    double t = 1.;
    //    mwIndex maxnnz = ceil(pow(log(n),3)*sqrt(n));
    //    maxnnz = std::min(maxnnz, (mwIndex)100000);
    mwIndex maxnnz = 10000;
    maxnnz = std::min(maxnnz, G.n);
 
    if (nrhs >= 5) { maxnnz = (mwIndex)mxGetScalar(prhs[4]);}
    if (nrhs >= 4) { t = mxGetScalar(prhs[3]); }
    if (nrhs >= 3) { eps = mxGetScalar(prhs[2]); }
    if ( maxnnz > G.n || maxnnz <= 0 ){
        mexErrMsgIdAndTxt("expmimv_mex:wrongParameterMaxnnz",
                          "expmimv_mex needs 1 <= maxnnz <= n ");
    }
    if ( t <= 0 ){
        mexErrMsgIdAndTxt("expmimv_mex:wrongParameterT",
                          "expmimv_mex needs 0 < t ");
    }

    if ( eps <= 0 || eps >= 1){
        mexErrMsgIdAndTxt("expmimv_mex:wrongParameterEps",
                          "expmimv_mex needs 0 < eps < 1 ");
    }
    
    std::vector< mwIndex > seeds;
    copy_array_to_index_vector( set, seeds );
    sparsevec hk;
    
    DEBUGPRINT(("expmimv_mex: preprocessing end: \n"));
    
    expm_svec(&G, seeds, hk, t, eps, maxnnz, npushes);
    
    DEBUGPRINT(("expmimv_mex: call to expm_svec() done\n"));
    
    if (nlhs > 0) { // sets output "hk" to the heat kernel vector computed
        mxArray* hkvec = mxCreateDoubleMatrix(G.n,1,mxREAL);
        plhs[0] = hkvec;
        double *ci = mxGetPr(hkvec);
        for (sparsevec::map_type::iterator it=hk.map.begin(),itend=hk.map.end();
             it!=itend;++it) {
            ci[it->first] = it->second;
        }
    }
}
