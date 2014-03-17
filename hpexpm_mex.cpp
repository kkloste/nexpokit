/**
 * @file hpexpm_mex.cpp
 * @author Kyle Kloster, David F. Gleich
 */

/**
 *
 * This file implements exp(P)*e_c via an approximate matrix-vector method
 * that uses a heap to determine the largest magnitude entries
 * of the iterative vector, then using only those nonzeros
 * to do the matvec.
 *
 * A - input matrix, must be sparse, n-by-n, and have nonnegative entries.
 *      This code is written specifically with column-stochastic A in mind,
 *      though it will work for non-stochastic matrices (probably poorly).
 *
 * The 'tol' input is intended for the case that the resulting
 * matvec's 1-norm is known a priori, so the function can quit adding
 * entries to the heap once it "weighs enough". (within 'tol'
 * of the 1-norm of the exactly computed matvec).
 *
 * TOL is temporarily removed -- I'll finish that aspect of the code later
 *  10/10/13, Kyle
 *
 * The 'maxnnz' input determines how many nonzero entries
 * of the input vector you want to use, at most.
 *
 */

#include "mex.h"
#include <vector>
#include <assert.h>
#include <math.h>
#include "minheapnoL.hpp" // include our heap functions


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

/**
 * @param n - sparse matrix size
 * @param cp - sparse matrix column pointer (CSC)
 * @param ari - sparse matrix row index (CSC)
 * @param a - sparse matrix values (CSC)
 * @param c - column index
 * @param maxnnz - the maximum number of nnz from input vec to use (1 <= maxnnz <= n)
 * @param y - the output vector (length n)
 */

void hpexpm(sparserow* G,
            const mwIndex c, const double eps, const double t, const mwIndex maxnnz,
            double* y)
{
    // Note: N is taken to be of type "double" because
    // it's used in computations below as type "double".
    
    mwIndex N = taylordegree(t,eps);
    mwIndex n = G->n;

    mwSize listsize = 0;
    mwSize maxlistsize = std::min( 4*maxnnz, n);
    // allocate data
    
    std::vector<mwIndex> Tvec(maxnnz,0);
    std::vector<double> Vvec(maxnnz,0);
    mwIndex *T = &Tvec[0];
    double *v = &Vvec[0];
	
    std::vector<mwIndex> list(maxlistsize,0);
    
    // mexPrintf("pass initialize \n");
    // Unroll first iterate of Horner's Rule
    for (mwIndex nzi=G->ai[c]; nzi < G->ai[c+1]; ++nzi) {
		mwIndex arinzi = G->aj[nzi];
        y[arinzi] += G->a[nzi]/N;
		list[listsize] = arinzi;
		listsize++;
    }
	if (y[c]==0){
		list[listsize] = c;
		listsize++;
	}
    y[c] = y[c] + 1;
    //	mexPrintf("pass first iteration \n");
    /***
     * STEP 1: BUILD HEAP FOR v
     *
     * STEP 2: Do heap-matvec y = A*v
     *
     ***/
    
    
    
    // HORNER'S RULE ITERATES
    for (mwIndex k=1; k <= N-1 ; k++){
		mwIndex Nk = (mwIndex)N-k;
        mwIndex listind = 0;
        mwIndex hsize = 0;

        // for each entry of y, check if we put it in the heap
        while ( listind < listsize ){
            mwIndex ind = list[listind];
            double valind = y[ind];
            
            // if the value is nonzero, we might add it to heap
            if (valind > 0) {
                // if the heap is full, we might have to delete
                if ( hsize >= maxnnz ){
                    double minval = y[T[0]];
                    // if new entry is big enough, add it
                    if ( minval < valind ){
                        // replace T[0] with y[ind], then heap_down
                        mwIndex ri = T[0];
                        T[0] = ind;
                        y[ri] = 0; // drop that entry from y
                        heap_down(0, hsize, T, y);
                        minval = y[T[0]];
                    }
                    else y[ind] = 0; // y[ind] so small it's excluded from heap
                }
                // if heap is not full, just place y[ind] at back of heap, then heap_up
                else {
                    T[hsize] = ind;
                    hsize++;
                    heap_up(hsize-1, hsize, T, y);
                }
            }
            listind ++;
        }
        // Copy the nonzeros in the heap of y into v
        // and zero out y, so it can be set equal to
        // the matvec y = A*v
        for (mwIndex ind = 0; ind < hsize ; ind++){
            mwIndex Tind = T[ind];
            v[ind] = y[Tind];
            y[Tind] = 0;
        }
        
        // At this point, y[] = 0 entirely, and v contains
        // just the entries from y that survived the heap process,
        // i.e. the maxnnz largest entries that were in y.
        
        listsize = 0;
        for (mwIndex ind = 0; ind < hsize; ind++) {
            mwIndex Tind = T[ind];
            double valind = v[ind]/(Nk);
            for ( mwIndex nzi=G->ai[Tind]; nzi < G->ai[Tind+1]; ++nzi) {
                mwIndex arinzi = G->aj[nzi];
                if (y[arinzi] == 0){
                    list[listsize] = arinzi;
                    listsize++;
                    if (listsize >= maxlistsize){ //resize list if it overflows
                        maxlistsize = maxlistsize*4;
                        list.resize(maxlistsize,0);
                    }
                }
                y[arinzi] += G->a[nzi]*valind;
            }
            v[ind] = 0;
        }
        if (y[c] == 0){
			list[listsize] = c;
			listsize++;
            if (listsize >= maxlistsize){ //resize list if it overflows
                maxlistsize = maxlistsize*4;
                list.resize(maxlistsize,0);
            }
		}
        y[c] = y[c] + 1;
		
        //		mexPrintf("iteration \n");
    }//terms of Taylor are complete
    return;
}

// USAGE  y = hpexpm_mex(A,c,eps,t,maxnnz)
void mexFunction(
                 int nargout, mxArray *pargout[],       // these are your outputs
                 int nargin, const mxArray *pargin[])   // these are your arguments
{
    // inputs
    // A - sparse n-by-n
    // c - positive integer, the column index of exp(A) to be computed
    // eps - accuracy desired
    // t - scalar, e^(tA)
    // maxnnz - integer scalar max number of entries to use in the matvecs
    
    const mxArray* mat = pargin[0];
    if ( mxIsSparse(mat) == false ){
        mexErrMsgIdAndTxt("hpexpm_mex:wrongInputMatrix",
                          "hpexpm_mex needs sparse input matrix");
    }
    if ( mxGetM(mat) != mxGetN(mat) ){
        mexErrMsgIdAndTxt("hpexpm_mex:wrongInputMatrixDimensions",
                          "hpexpm_mex needs square input matrix");
    }

    // decode the sparse matrix
    sparserow G;
    G.m = mxGetM(mat);
    G.n = mxGetN(mat);
    G.ai = mxGetJc(mat);
    G.aj = mxGetIr(mat);
    G.a = mxGetPr(mat);

    mwIndex c = (mwIndex)mxGetScalar(pargin[1])-1;
    
    double eps = 1e-4;
    double t = 1.;
    mwIndex maxnnz = std::min((mwIndex)1000,G.n);
    if (nargin >= 3){ eps = mxGetScalar(pargin[2]); }
    if (nargin >= 4){ t = mxGetScalar(pargin[3]); }
    if (nargin >= 5){ maxnnz = (mwIndex)mxGetScalar(pargin[4]); }
    
	
    pargout[0] = mxCreateDoubleMatrix(G.n,1,mxREAL);
    
    double* y = mxGetPr(pargout[0]);
    
    // mexPrintf("Starting call \n");
    
    if ( nargin > 5 || nargin < 2 ){
        mexErrMsgIdAndTxt("hpexpm_mex:wrongNumberInputs",
                          "hpexpm_mex needs 3 inputs");
    }
    if ( c < 0 || c >= G.n ){
        mexErrMsgIdAndTxt("hpexpm_mex:wrongParamterC",
                          "hpexpm_mex needs 1 <= c <= n");
    }
    if ( maxnnz > G.n || maxnnz < 1 ){
        mexErrMsgIdAndTxt("hpexpm_mex:wrongParamterMaxnnz",
                          "hpexpm_mex needs 1 <= maxnnz <= n");
    }
    
    
    hpexpm(&G, // sparse matrix
           c, eps, t, maxnnz, // parameters
           y); // output product vector
}