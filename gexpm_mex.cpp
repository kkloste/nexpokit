/** 
 * @file gexpm_mex.cpp
 * @author Kyle Kloster, David F. Gleich
 */

/**
 * This file implements a gauss-southwell type method for the truncated 
 * taylor series approximation for a column of the matrix exponential
 */

#include "mex.h"

#include "heap.hpp" // include our heap functions

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
 */
void gexpm(const mwSize n, const mwIndex* cp, const mwIndex* ari, const double* a, 
           const mwIndex c, const mwIndex N,  const double tol, const mwIndex maxsteps, 
           double* y)
{
    //mexPrintf("Input n=%i N=%i c=%i tol=%i maxsteps=%i\n", n, N, c, tol, maxsteps);
    mwIndex M = n*N;
    double sumresid = 0.;
    
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
    r[rentry(c,0)] = 1;
    sumresid += 1.;
    T[hsize] = rentry(c,0);
    L[rentry(c,0)] = hsize;
    hsize++;
    heap_up(hsize-1, hsize, T, L, r);
    
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
        mwIndex ri = T[0];
        T[0] = T[hsize-1];
        L[T[0]] = 0;
        L[ri] = M+1; /* the null location, the order here is important! */
        hsize--;
        // mexPrintf("step %05i - popped ri=%3i i=%3i j=%3i rij=%.18e\n",
        //             iter, ri, ri%n, ri/n, r[ri]);
        heap_down(0, hsize, T, L, r);
        
        // decode incides i,j
        mwIndex i = ri%n;
        mwIndex j = ri/n;
        
        double rij = r[ri];
        
        // update yi
        y[i] += rij;
        
        // update r, no need to update heap here 
        r[ri] = 0;
        sumresid -= rij;
        double rijs = rij/(double)(j+1);
        
        if (j == N-1) {
            // this is the terminal case, and so we add the column of A 
            // directly to the solution vector y
            for (mwIndex nzi=cp[i]; nzi < cp[i+1]; ++nzi) {
                mwIndex v = ari[nzi];
                double ajv = a[nzi];
                y[v] += ajv*rijs;
            }
        } else {
            // this is the interior case, and so we add the column of A 
            // to the residual at the next time step.
            for (mwIndex nzi=cp[i]; nzi < cp[i+1]; ++nzi) {
                mwIndex v = ari[nzi];
                double ajv = a[nzi];
                mwIndex re = rentry(v,j+1);
                r[re] += ajv*rijs;
                sumresid += ajv*rijs;
                
                if (L[rentry(v,j+1)] > M) { // then this is a new heap entry
                    T[hsize] = re;
                    L[re] = hsize;
                    hsize ++;
                }
                mwIndex k = L[re];
                k = heap_down(k, hsize, T, L, r);
                heap_up(k, hsize, T, L, r);
            }
        }
        
        
        if (sumresid < tol) {
            break;
        }
    }
    return; // because we "break" out of for loop
}

void mexFunction(
  int nargout, mxArray *pargout[],       // these are your outputs
  int nargin, const mxArray *pargin[])   // these are your arguments
{
    // inputs
    // A - sparse n-by-n
    // c - integer scalar 1 <= c <= n
    // N - integer scalar 1 <= N <= Inf
    // tol - double value, 0 < tol < Inf
    // maxsteps - integer scalar max-steps 
    
    const mxArray* A = pargin[0];
    mwIndex c = (mwIndex)mxGetScalar(pargin[1])-1;
    mwIndex N = (mwIndex)mxGetScalar(pargin[2]);
    double tol = mxGetScalar(pargin[3]);
    mwIndex maxsteps = (mwIndex)mxGetScalar(pargin[4]);
    
    mwSize n = mxGetM(A);
    pargout[0] = mxCreateDoubleMatrix(n,1,mxREAL);
    
    // decode the sparse matrix
    mwIndex* cp = mxGetJc(A);
    mwIndex* ri = mxGetIr(A);
    double* a = mxGetPr(A);
    
    double* y = mxGetPr(pargout[0]);
    
    // mexPrintf("Starting call \n");
    
    mxAssert(N > 0, "N must be bigger than 0");
    mxAssert(tol > 0 && tol <= 1, "tol must be 0 < tol <= 1");
    mxAssert(c >= 0 && c < n, "column c must be 1 <= c <= n");
    mxAssert(maxsteps > 0, "we must have maxsteps >= 0");
    
    gexpm(n, cp, ri, a, // sparse matrix
          c, N, tol, maxsteps, // parameters
          y);
    
}
    
    
  