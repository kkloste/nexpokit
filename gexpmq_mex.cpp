/** 
 * @file gexpmq_mex.cpp
 * @author Kyle Kloster, David F. Gleich
 */

/**
 *  This improves on gexpmq_mex by scaling block-k of the residual
 *  by psi_k(t) for all k -- gexpmq simply ran until every entry
 *  of the residual  satisfied r(i,j) < tol/(n*N). Instead, this
 *  runs until entries in r(:,j) are < tol/(n*psivec[j]), which
 *  saves some computation and is more easily adapted to e^(tP)
 *  (whereas t=1 must be the case for gexpmq)
 */

#ifndef __APPLE__
#define __STDC_UTF_16__ 1
#endif 

#include <vector>
#include <queue>
#include <utility> // for pair sorting
#include <assert.h>
#include <limits>
#include <algorithm>
#include <math.h>
#include <sparsehash/dense_hash_map>

#include "taydeg.hpp"

#include <mex.h>

#define DEBUGPRINT(x) do { if (debugflag) { \
mexPrintf x; mexEvalString("drawnow"); } \
} while (0)

int debugflag = 0;


#include <sys/types.h>
#include <sys/timeb.h>
#include <sys/time.h>
double sf_time()
{
#if defined(_WIN32) || defined(_WIN64)
  struct __timeb64 t; _ftime64(&t);
  return (t.time*1.0 + t.millitm/1000.0);
#else
  struct timeval t; gettimeofday(&t, 0);
  return (t.tv_sec*1.0 + t.tv_usec/1000000.0);
#endif
}

struct sparsevec {
    typedef google::dense_hash_map<mwIndex,double> map_type;
    map_type map;

    sparsevec() {
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

struct local_stochastic_graph_exponential
{
    // inputs
    sparserow* G;
    double t, eps;
    size_t maxpush; // used to cap work below convergence
    
    // derived
    mwIndex n, N;
    std::vector<double> pushcoeff;
    
    // the local graphindex
    typedef unsigned int lindex; 
    
    // sentinal
    const lindex sval; 
    
    // the local graph
    google::dense_hash_map<mwIndex, lindex > imap;
    lindex nextind;
    lindex noffset; // offset is next index of lneigh
                    // that is unused
    
    std::vector< lindex > lneigh;                  // local set of neigbhors
    std::vector< std::pair< lindex, lindex > > lp; // local set of offsets
    std::vector< mwIndex > gmap;                   // global indices
    
    // outputs
    size_t npush;
    
    // data
    // chosen so initial data is about 256k
    const static int nstart = 2048;

      
    local_stochastic_graph_exponential(
            sparserow* G_, 
            const double t_, double eps_,
            size_t maxpush_) 
	: G(G_), t(t_), eps(eps_), maxpush(maxpush_), 
            n(G->n), N(taylordegree(t, eps)),
            pushcoeff(N+1,0.), sval(std::numeric_limits<lindex>::max()),
            nextind(0), noffset(0),
            lneigh(nstart, sval), 
            lp(nstart, std::make_pair(sval, sval)), 
            gmap(nstart, (mwIndex)-1), npush(0)
    {        
        DEBUGPRINT(("gsqres interior: t=%lf eps=%lf \n", t, eps));
        DEBUGPRINT(("gsqres: n=%i N=%i \n", n, N));
        
        // initialize the weights for the different residual partitions
        //  into the vector "pushcoeff"
        std::vector<double> psivec(N+1,0.);
        psivec[N] = 1;
        for (lindex k = 1; k <= N ; k++){
            psivec[N-k] = psivec[N-k+1]*t/(double)(N-k+1) + 1;
        } // psivec[k] = psi_k(t)

        pushcoeff[0] = ((exp(t)*eps)/(double)N)/psivec[0];
        // This is a more numerically stable way to compute
        //      pushcoeff[j] = exp(t)*eps/(n*N*psivec[j])
        for (lindex k = 1; k <= N ; k++){
            pushcoeff[k] = pushcoeff[k-1]*(psivec[k-1]/psivec[k]);
        }
        
        if (n > std::numeric_limits<lindex>::max()) {
            mexErrMsgIdAndTxt("gexpmq_mex:unimplemented", 
                    "support only up to %i elements", 
                    std::numeric_limits<lindex>::max());
        }
        
        imap.set_empty_key(sval);        
    }
    
    template <class T>
    void grow_vector_to_index(std::vector< T >& vec, lindex ind, T val) {
        if (ind >= vec.size()) {
            vec.resize(ind+1, val); 
        }
    }
    
    void fetch_row(lindex li) {
        if (li < lp.size() && lp[li].first != sval) {
            return;
        }
        
        mwIndex gi = gmap[li];
        
        // make sure we have enough room for this one
        grow_vector_to_index(lp, li, std::make_pair(sval,sval));
        
        lindex sindex = noffset;
        lindex eindex = noffset;
        lindex deg = sr_degree(G,gi);
        grow_vector_to_index(lneigh, sindex+deg-1, sval);
        
        // copy everything to the local graph
        for (mwIndex nzi=G->ai[gi]; nzi < G->ai[gi+1]; ++nzi) {
            lneigh[eindex] = fetch_index(G->aj[nzi]);
            eindex ++;
        }
        lp[li] = std::make_pair(sindex, eindex);
        noffset = eindex;
    }
    
    // fetch the row associated with the given index 
    // and assign it if it doesn't exist
    lindex fetch_index(mwIndex gi) {
        if (imap.count(gi) != 0) {
            return imap[gi];
        }
        lindex li = nextind;
        nextind ++;
        imap[gi] = li;
        grow_vector_to_index(gmap, li, (mwIndex)(-1));
        gmap[li] = gi;
        return li;
    }
    

    mwIndex compute(std::vector<mwIndex>& set, sparsevec& yout) {
        
        // allocate residuals for each level
        std::vector< std::vector< double > > resid(N);
        std::vector< double > y;

        // this is defined earlier, but outside the scope of compute(), so I redefine it here
        std::vector<double> psivec(N+1,0.);
        psivec[N] = 1;
        for (lindex k = 1; k <= N ; k++){
            psivec[N-k] = psivec[N-k+1]*t/(double)(N-k+1) + 1;
        } // psivec[k] = psi_k(t)
     
        double sumresid = 0.;
        // the queue stores the entries for the next time step
        std::queue< lindex > Q;
        
        // set the initial residual, add to the queue
        for (size_t i=0; i<set.size(); ++i) {
            mwIndex ri = set[i];
            lindex li = fetch_index(ri);
            //rij = value in entry ri of the input vector
            double rij = 1.;
            sumresid += rij*psivec[0];;
            grow_vector_to_index(resid[0], li, 0.);
            resid[0][li] += rij;
            Q.push( li );
        }
        
        for (lindex j = 0; j < N; ++j) {
            // STEP 0: determine how many entries are in the queue
            // and determine what the residual tolerance is to 
            // push
            size_t qsize = Q.size();
            double pushtol = pushcoeff[j]/(double)qsize;
            
            // STEP 1: process each of the next elements in the 
            // queue and add them to the next level
            
            for (size_t qi = 0; qi < qsize; ++qi) {
                lindex i = Q.front();
                Q.pop();
                            
                double rij = resid[j][i];
            
                if (rij < pushtol) {
                    // skip to the next iteration of the loop
                    continue;
                }
            
                // STEP 2: fetch the row
                fetch_row(i);
            
                // find the degree
                std::pair<lindex, lindex> offsets = lp[i];
                double degofi = (double)(offsets.second - offsets.first);
            
                // update the solution
                grow_vector_to_index(y, i, 0.);
                y[i] += rij;
            
                resid[j][i] = 0.;
                sumresid -= rij*psivec[j];;
            
                double rijs = t*rij/(double)(j+1);
                double ajv = 1./degofi;
                double update = rijs*ajv;
            
                if (j == N-1) {
                    // this is the terminal case, and so we add the column of A
                    // directly to the solution vector y
                    for (lindex nzi = offsets.first; nzi < offsets.second; ++nzi) {
                        lindex v = lneigh[nzi];
                        grow_vector_to_index(y, v, 0.);
                        y[v] += update;
                    }
                    npush += degofi;
                } else {
                    // this is the interior case, and so we add the column of A
                    // to the residual at the next time step.
                    for (lindex nzi = offsets.first; nzi < offsets.second; ++nzi) {
                        lindex v = lneigh[nzi];
                        grow_vector_to_index(resid[j+1], v, 0.);
                        double reold = resid[j+1][v]; 
                        double renew = reold + update;
                        sumresid += update*psivec[j+1];;
                        resid[j+1][v] = renew;
                        // For this implementation, we need all 
                        // non-zero residuals in the queue.
                        if (reold == 0.) {
                            Q.push( v );
                        }
                    }
                    npush += degofi;
                }
                
                if (maxpush > 0 && npush >= maxpush) { break; }
                if (sumresid < eps*exp(t)) { break; }
                // terminate when Q is empty, i.e. we've pushed all r(i,j) > exp(t)*eps*/(N*n*psi_j(t))
            }
            if (sumresid < eps*exp(t)) { break; }
            if (maxpush > 0 && npush >= maxpush) { break; }
        } // end 'for'
        
        
        // store the output back in the solution vector y
        for (lindex i=0; i<y.size(); ++i) {
            if (y[i] > 0) {
                yout.map[gmap[i]] = y[i];
            }
        }
        
        return npush;
    }  
};


    


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
// [y npushes] = gexpmq_mex(A,set,eps,t,debugflag)
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    if (nrhs < 2 || nrhs > 6) {
        mexErrMsgIdAndTxt("gexpmq_mex:wrongNumberArguments",
                          "gexpmq_mex needs two to five arguments, not %i", nrhs);
    }
    if ( nlhs > 2 ){
        mexErrMsgIdAndTxt("gexpmq_mex:wrongNumberOutputs",
                          "gexpmq_mex needs two outputs, not %i", nlhs);
    }
    
    if (nrhs == 5) {
        debugflag = (int)mxGetScalar(prhs[4]);
    }
    
    DEBUGPRINT(("gexpmq_mex: preprocessing start: \n"));
    
    const mxArray* mat = prhs[0];
    const mxArray* set = prhs[1];
    
    if ( mxIsSparse(mat) == false ){
        mexErrMsgIdAndTxt("gexpmq_mex:wrongInputMatrix",
                          "gexpmq_mex needs sparse input matrix");
    }
    if ( mxGetM(mat) != mxGetN(mat) ){
        mexErrMsgIdAndTxt("gexpmq_mex:wrongInputMatrixDimensions",
                          "gexpmq_mex needs square input matrix");
    }
    
    double npushesmem = 0.;
    double* npushes = &npushesmem;
    if (nlhs > 1){
        plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
        npushes = mxGetPr(plhs[1]);
    }
    
    double eps = pow(10,-5);
    double t = 1.;
    mwIndex maxpush = 0;

    if (nrhs >= 4) { t = mxGetScalar(prhs[3]); }
    if (nrhs >= 3) { eps = mxGetScalar(prhs[2]); }
    if (nrhs >= 6) { maxpush = (mwIndex) mxGetScalar(prhs[5]); }
    
    
    sparserow G;
    G.m = mxGetM(mat);
    G.n = mxGetN(mat);
    G.ai = mxGetJc(mat);
    G.aj = mxGetIr(mat);
    G.a = mxGetPr(mat);
    
    std::vector< mwIndex > seeds;
    copy_array_to_index_vector( set, seeds );
    sparsevec hk;
    
    DEBUGPRINT(("gexpmq_mex: preprocessing end: \n"));
    
    //gexpmq(&G, seeds, hk, t, eps, npushes);
    local_stochastic_graph_exponential c(&G, t, eps, (size_t)maxpush);
    *npushes = c.compute(seeds, hk);

    DEBUGPRINT(("gexpmq_mex: call to gsqres() done\n"));

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
