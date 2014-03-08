/** 
 * @file gsqres_mex.cpp
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

#include <vector>
#include <queue>
#include <utility> // for pair sorting
#include <assert.h>
#include <limits>
#include <algorithm>
#include <math.h>

#ifdef __APPLE__
#include <tr1/unordered_set>
#include <tr1/unordered_map>
#define tr1ns std::tr1
#else
#include <unordered_set>
#include <unordered_map>
#define __STDC_UTF_16__ 1
#define tr1ns std
#endif

#include <mex.h>

#define DEBUGPRINT(x) do { if (debugflag) { \
mexPrintf x; mexEvalString("drawnow"); } \
} while (0)

int debugflag = 0;

struct sparsevec {
    typedef tr1ns::unordered_map<mwIndex,double> map_type;
    map_type map;
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


/**
 * Computes the degree N for the Taylor polynomial
 * of exp(tP) to have error less than eps*exp(t)
 *
 * ( so exp(-t(I-P)) has error less than eps )
 */
unsigned int taylordegree(const double t, const double eps) {
    double eps_exp_t = eps*exp(t);
    double error = exp(t)-1;
    double last = 1.;
    double k = 0.;
    while(error > eps_exp_t){
        k = k + 1.;
        last = (last*t)/k;
        error = error - last;
    }
    return std::max((int)k, (int)1);
}

/*****
 *
 *          above:  DATA STRUCTURES
 *
 *
 *
 *          below:  EXPM FUNCTION
 *
 ****/

/**
 *
 *  gsqexpmseed inputs:
 *      G   -   adjacency matrix of an undirected graph
 *      set -   seed vector: the indices of a seed set of vertices
 *              around which cluster forms; normalized so
 *                  set[i] = 1/set.size(); )
 *  output:
 *      y = exp(tP) * set
 *              with relative 1-norm error of eps
 *  parameters:
 *      t   - the value of t
 *      eps - the accuracy
 */
int gexpmq(sparserow* G, std::vector<mwIndex>& set, sparsevec& y,
                const double t, const double eps)
{
    DEBUGPRINT(("gsqexpmseed interior: t=%f eps=%f \n", t, eps));
    mwIndex n = G->n;
    mwIndex N = (mwIndex)taylordegree(t, eps);
    DEBUGPRINT(("gsqexpmedseed: n=%i N=%i \n", n, N));
    
    // initialize the weights for the different residual partitions
    // r(i,j) > exp(t)*eps/(N*psi_j(t))
    //  into the vector "pushcoeff"
    std::vector<double> psivec(N+1,0.);
    psivec[N] = 1;
    for (int k = 1; k <= N ; k++){
        psivec[N-k] = psivec[N-k+1]*t/(double)(N-k+1) + 1;
    } // psivec[k] = psi_k(t)
    std::vector<double> pushcoeff(N+1,0.);
    pushcoeff[0] = (((exp(t)/(double)n)*eps)/(double)N)/psivec[0];
    // This is a more numerically stable way to compute
    //      pushcoeff[j] = exp(t)*eps/(n*N*psivec[j])
    for (int k = 1; k <= N ; k++){
        pushcoeff[k] = pushcoeff[k-1]*(psivec[k-1]/psivec[k]);
    }
    
    mwIndex ri = 0;
    mwIndex npush = 0;
    double rij = 0;
    // allocate data
    sparsevec rvec;
    
    // i is the node index, j is the "step"
#define rentry(i,j) ((i)+(j)*n)
    std::queue<mwIndex> Q;
    // set the initial residual, add to the queue
    for (size_t i=0; i<set.size(); ++i) {
        ri = set[i];
        //rij = value in entry ri of the input vector
        rij = 1.;
        rvec.map[rentry(ri,0)]+=rij;
        Q.push(rentry(ri,0));
    }
    
//    while (npush < max_push_count) {
    while (1) {
        // STEP 1: pop top element off of heap
        ri = Q.front();
        Q.pop();
        // decode incides i,j
        mwIndex i = ri%n;
        mwIndex j = ri/n;
        
        double degofi = (double)sr_degree(G,i);
        rij = rvec.map[ri];
        
        // update yi
        y.map[i] += rij;
        
        // update r, no need to update heap here
        rvec.map[ri] = 0;
        
        double rijs = t*rij/(double)(j+1);
        double ajv = 1./degofi;
        double update = rijs*ajv;
        
        if (j == N-1) {
            // this is the terminal case, and so we add the column of A
            // directly to the solution vector y
            for (mwIndex nzi=G->ai[i]; nzi < G->ai[i+1]; ++nzi) {
                mwIndex v = G->aj[nzi];
                y.map[v] += update;
            }
            npush += degofi;
        }
        else {
            // this is the interior case, and so we add the column of A
            // to the residual at the next time step.
            for (mwIndex nzi=G->ai[i]; nzi < G->ai[i+1]; ++nzi) {
                mwIndex v = G->aj[nzi];
                mwIndex re = rentry(v,j+1);
                double reold = rvec.get(re);
                double renew = reold + update;
                rvec.map[re] = renew;
                if (renew >= pushcoeff[j+1] && reold < pushcoeff[j+1]) {
                    Q.push(re);
                }
            }
            npush+=degofi;
        }
        // terminate when Q is empty, i.e. we've pushed all r(i,j) > exp(t)*eps*/(N*n*psi_j(t))
        if ( Q.size() == 0) { return npush; }
    }//end 'while'
    return (npush);
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
// [y npushes] = gsqres_mex(A,set,eps,t,debugflag)
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    if (nrhs < 2 || nrhs > 5) {
        mexErrMsgIdAndTxt("gexpmq_mex:wrongNumberArguments",
                          "gexpmq_mex needs two to five arguments, not %i", nrhs);
    }
    mxAssert(nlhs <= 2, "Too many output arguments");
    if (nrhs == 5) {
        debugflag = (int)mxGetScalar(prhs[4]);
    }
    DEBUGPRINT(("gexpmq_mex: preprocessing start: \n"));
    
    const mxArray* mat = prhs[0];
    const mxArray* set = prhs[1];
    
    mxAssert(mxIsSparse(mat), "Input matrix is not sparse");
    mxAssert(mxGetM(mat) == mxGetN(mat), "Input matrix not square");
    mxArray* npushes = mxCreateDoubleMatrix(1,1,mxREAL);
    
    double eps = pow(10,-5);
    double t = 1.;

    if (nrhs >= 4) { t = mxGetScalar(prhs[3]); }
    if (nrhs >= 3) { eps = mxGetScalar(prhs[2]); }
    
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
    
    *mxGetPr(npushes) = (double)gexpmq(&G, seeds, hk, t, eps);

    DEBUGPRINT(("gexpmq_mex: call to gexpmq() done\n"));

    if (nlhs > 1) { plhs[1] = npushes; }
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