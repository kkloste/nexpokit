/** 
 * @file gexpm_mex.cpp
 * @author Kyle Kloster, David F. Gleich
 */

/**
 * This file implements a gauss-southwell type method for the truncated 
 * taylor series approximation for a column of the matrix exponential
 */
#ifndef __APPLE__
#define __STDC_UTF_16__ 1
#endif 

#include "mex.h"
#include <math.h>
#include "heap.hpp" // include our heap functions
#include "taydeg.hpp" // for computing Taylor degree

#include <vector>


/**
 * @param n - sparse matrix size
 * @param cp - sparse matrix column pointer (CSC)
 * @param ari - sparse matrix row index (CSC)
 * @param a - sparse matrix values (CSC)
 * @param c - the column index to approximate (0 <= c < n)
 * @param N - the number of steps (1 <= N < Inf)
 * @param tol - the stopping tolerance (0 < tol < Inf)
 * @param maxsteps - the maximum number of steps to take (1 <= maxsteps < Inf)
 * @param y - the output vector (length n)
 * @param nsteps - the number of output steps (length 1)
 */
void gexpm(const mwSize n, const mwIndex* cp, const mwIndex* ari, const double* a, 
           const mwIndex c, const mwIndex N,  const double tol, const mwIndex maxsteps, 
           double* y, double *nsteps, double *npushes)
{
    //mexPrintf("Input n=%i N=%i c=%i tol=%i maxsteps=%i\n", n, N, c, tol, maxsteps);

    std::vector<double> psivec(N+1,0.);
    psivec[N] = 1;
    double t = 1.;
    for (mwIndex k = 1; k <= N ; k++){
        psivec[N-k] = psivec[N-k+1]*t/(double)(N-k+1) + 1;
    } // psivec[k] = psi_k(t)
    
    mwIndex M = n*N;
    // allocate data
    std::vector<double> rvec(M,0.);
    double *r = &rvec[0];
    
    mwIndex hsize = 0;
    std::vector<mwIndex> Lvec(M,(M+1));
    std::vector<mwIndex> Tvec(M,0);
    mwIndex *L = &Lvec[0];
    mwIndex *T = &Tvec[0];
    
    // i is the node index, j is the "step"
    #define rentry(i,j) ((i)+(j)*n)
    
    // mexPrintf("Init...\n");
 
    // set the initial residual, add to the heap, and update
    double sumresid = 0.;
    r[rentry(c,0)] = 1;
    sumresid += 1.;
    T[hsize] = rentry(c,0);
    L[rentry(c,0)] = hsize;
    hsize++;
    heap_up(hsize-1, hsize, T, L, r);
    
    mwIndex npush = 0;
    *nsteps = (double)maxsteps; // set the default, which we change on early exit
    
    mwIndex i,j,v,re,k,ri;
    double rijs, rij, ajv;
    double toln = tol/(double)n;
    *nsteps = (double)maxsteps; // set the default, which we change on early exit
    
    //mexPrintf("Loop...\n");
    
    for (mwIndex iter = 0; iter < maxsteps; ++iter) {
        
        /* STEP 1: pop top element off of heap
         *  * get indices i,j from T
         *  * add r(i,j) to y(i)
         *  * set r(i,j) to zero (update sumresid)
         * STEP 2: get i^th column from A
         *  * get neighbors of ith node
         *  * (if j == N-1), add the column to y instead of r.
         *  * add as a column to next time-step of r, and update heap
         *  *  (update sumresid)
         * Check for convegence!
        */
        
        // STEP 1: pop top element off of heap
        ri = T[0];
        T[0] = T[hsize-1];
        L[T[0]] = 0;
        L[ri] = M+1; /* the null location, the order here is important! */
        hsize--;
        // mexPrintf("step %05i - popped ri=%3i i=%3i j=%3i rij=%.18e\n",
        //             iter, ri, ri%n, ri/n, r[ri]);
        heap_down(0, hsize, T, L, r);
        
        // decode incides i,j
        i = ri%n;
        j = ri/n;
        
        rij = r[ri];
        
        // update yi
        y[i] += rij;
        
        // update r, no need to update heap here 
        r[ri] = 0;
        sumresid -= rij;
        rijs = rij/(double)(j+1);
        
        if (j == N-1) {
            // this is the terminal case, and so we add the column of A 
            // directly to the solution vector y
            for (mwIndex nzi=cp[i]; nzi < cp[i+1]; ++nzi) {
                v = ari[nzi];
                ajv = a[nzi];
                y[v] += ajv*rijs;
            }
        } else {
            // this is the interior case, and so we add the column of A 
            // to the residual at the next time step.
            for (mwIndex nzi=cp[i]; nzi < cp[i+1]; ++nzi) {
                v = ari[nzi];
                ajv = a[nzi];
                re = rentry(v,j+1);
                r[re] += ajv*rijs;
                sumresid += ajv*rijs;
                
                if (r[re] > toln/psivec[j+1]) {
                    if (L[rentry(v,j+1)] > M) { // then this is a new heap entry
                        T[hsize] = re;
                        L[re] = hsize;
                        hsize ++;
                    }
                    k = L[re];
                    k = heap_down(k, hsize, T, L, r);
                    heap_up(k, hsize, T, L, r);
                }
            }
            npush+=cp[i+1]-cp[i];
        }
        if (sumresid < tol) {
            *nsteps = (double)iter;
            break;
        }
    }
    
    *npushes = (double)npush;
    
    return; // because we "break" out of for loop
}

void mexFunction(
  int nargout, mxArray *pargout[],       // these are your outputs
  int nargin, const mxArray *pargin[])   // these are your arguments
{
    // inputs
    // A - sparse n-by-n
    // c - integer scalar 1 <= c <= n
    // tol - double value, 0 < tol < Inf
    // N - integer scalar 1 <= N <= Inf
    // maxsteps - integer scalar max-steps 
    
    const mxArray* A = pargin[0];
    mwIndex c = (mwIndex)mxGetScalar(pargin[1])-1;
    double tol = mxGetScalar(pargin[2]);
    mwSize n = mxGetM(A);
    mwIndex maxsteps = 10*n;
    if (nargin >= 4){
        maxsteps = (mwIndex)mxGetScalar(pargin[3]);
    }
    mwIndex N = (mwIndex)taylordegree(1.,tol);
    if (nargin == 5){
        N = (mwIndex)mxGetScalar(pargin[4]);
    }

    
    pargout[0] = mxCreateDoubleMatrix(n,1,mxREAL);
    pargout[1] = mxCreateDoubleMatrix(1,1,mxREAL);
    pargout[2] = mxCreateDoubleMatrix(1,1,mxREAL);
    
    // decode the sparse matrix

    mwIndex* cp = mxGetJc(A);
    mwIndex* ri = mxGetIr(A);
    double* a = mxGetPr(A);
    
    double* y = mxGetPr(pargout[0]);
    double* nsteps = mxGetPr(pargout[1]);
    double* npushes = mxGetPr(pargout[2]);
    
    // mexPrintf("Starting call \n");
    
    if ( N <= 0 ){
        mexErrMsgIdAndTxt("gexpm_mex:wrongInputParamterN",
                          "gexpm_mex needs N > 0");
    }
    if ( tol <= 0 || tol > 1){
        mexErrMsgIdAndTxt("gexpm_mex:wrongInputParamterTol",
                          "gexpm_mex needs 0 < tol <= 1");
    }
    if ( c >= n ){
        mexErrMsgIdAndTxt("gexpm_mex:wrongInputParamterC",
                          "gexpm_mex needs 1 <= c <= n");
    }
    if ( maxsteps <= 0 ){
        mexErrMsgIdAndTxt("gexpm_mex:wrongInputParamterN",
                          "gexpm_mex needs maxsteps > 0");
    }
    
    gexpm(n, cp, ri, a, // sparse matrix
          c, N, tol, maxsteps, // parameters
          y, nsteps, npushes);
    
}
    
    
  
