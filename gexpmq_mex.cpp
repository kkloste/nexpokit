/** 
 * @file gexpmq_mex.cpp
 * @author Kyle Kloster, David F. Gleich
 */

/**
 * This file implements an approximate gauss-southwell method using a Queue
 * instead of a heap to approximate the largest element for the truncated 
 * taylor series approximation for a column of the matrix exponential
 */
#ifndef __APPLE__
#define __STDC_UTF_16__ 1
#endif 

#include <queue>
#include <vector>
#include <utility> // for pair sorting
#include <assert.h>
#include <math.h>
#include "taydeg.hpp"

#include "mex.h"

#define DEBUGPRINT(x) do { if (debugflag) { \
mexPrintf x; mexEvalString("drawnow"); } \
} while (0)

int debugflag = 0;

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
void gexpmq(const mwSize n, const mwIndex* cp, const mwIndex* ari, const double* a, 
            const mwIndex c, const mwIndex N,  const double tol, const mwIndex maxsteps, 
            double* y, double *nsteps, double *npushes)
{
DEBUGPRINT(("Input n=%i N=%i c=%i tol=%i maxsteps=%i\n", n, N, c, tol, maxsteps));
    mwIndex M = n*N;
    double sumresid = 0.;
    double tolnN = tol/(n*N);

	// compute the tolerances for the components of the residual
    std::vector<double> psivec(N+1,0.);
    psivec[N] = 1;
    for (mwIndex k = 1; k <= N ; k++){
        psivec[N-k] = psivec[N-k+1]*t/(double)(N-k+1) + 1;
    } // psivec[k] = psi_k(t)

    pushcoeff[0] = (exp(t)*eps/(double)N)/psivec[0];
	// This is a more numerically stable way to compute
	//      pushcoeff[j] = exp(t)*eps/(N*psivec[j])
	for (lindex k = 1; k <= N ; k++){
		pushcoeff[k] = pushcoeff[k-1]*(psivec[k-1]/psivec[k]);
	}
	


    // allocate data
    
    std::vector<double> rvec(M,0.);
    double *r = &rvec[0];

DEBUGPRINT(("delcare queue\n"));
    std::queue<mwIndex> Q;
    
    // i is the node index, j is the "step"
    #define rentry(i,j) ((i)+(j)*n)
    
DEBUGPRINT(("Init...\n"));
 
    // set the initial residual, push onto queue
    r[rentry(c,0)] = 1;
    sumresid += 1.;
    Q.push(rentry(c,0));
    
    mwIndex npush = 0;
    *nsteps = (double)maxsteps; // set the default, which we change on early exit
    mwIndex iter = 0;
DEBUGPRINT(("Before for-loop\n"));
	for (mwIndex j = 0; j <= N ; j++){
		mwIndex numnnz = Q.size();
		double pushtol = pushcoeff[j]/(double)numnnz;
		for ( mwIndex Qentries = 1; <= numnnz ; Qentries++){
			if (iter > maxnnz){
				j = N+1;
				break;
			}
			iter++;
		
			// STEP 1: pop top element off of residual
			mwIndex ri = Q.front();
			Q.pop();
			double rij = r[ri];
			r[ri] = 0;
			sumresid -= rij;
			  
			// decode incides i,j
			mwIndex i = ri%n;
//			mwIndex j = ri/n;

			// update yi
			y[i] += rij;

			if ( rij >= pushtol ){ // not small enough to round it off
				double rijs = rij/(double)(j+1);
				if (j == N-1) {
					// this is the terminal case, and so we add the column of A 
					// directly to the solution vector y
					for (mwIndex nzi=cp[i]; nzi < cp[i+1]; ++nzi) {
						mwIndex v = ari[nzi];
						double ajv = a[nzi];
						y[v] += ajv*rijs;
					}
					npush+=cp[i+1]-cp[i];
				} else {
					// this is the interior case; add the column of A 
					// to the residual at the next time step.
					for (mwIndex nzi=cp[i]; nzi < cp[i+1]; ++nzi) {
						mwIndex v = ari[nzi];
						double ajv = a[nzi];
						mwIndex re = rentry(v,j+1);
						double reold = r[re];
						r[re] += ajv*rijs;
						sumresid += ajv*rijs;
						// if this wasn't in the queue, push it in
						if (r[re] > 0 && reold == 0) { Q.push(re); }
					}
					npush+=cp[i+1]-cp[i];
				}
			}
			if (sumresid < tol || Q.size() == 0) {
				Qentries = numnnz+1;
				j = N+1;
				break;
			}			
		} //end pushes for particular Taylor term
    } // end Taylor terms
	*nsteps = (double)iter;    
    *npushes = (double)npush;
    return; // because we "break" out of for loop
}

// USAGE: [expmv nstep npush] = gexpmq_mex(A,c,tol)
// optional: gexpmq_mex(A,c,tol,maxsteps,N)
void mexFunction(
  int nargout, mxArray *pargout[],       // these are your outputs
  int nargin, const mxArray *pargin[])   // these are your arguments
{
    // inputs
    // A - sparse n-by-n
    // c - integer scalar 1 <= c <= n
    // tol - double value, 0 < tol < Inf
    // maxsteps - integer scalar max-steps 
    // optional: N - integer scalar 1 <= N <= Inf
    // if N is not given, function 'tayordegree' is
    // called to set N
    const mxArray* A = pargin[0];
    mwSize n = mxGetM(A);
    mwIndex c = (mwIndex)mxGetScalar(pargin[1])-1;
    double tol = mxGetScalar(pargin[2]);
    mwIndex maxsteps = 10*n;
    if (nargin >= 4){
        maxsteps = (mwIndex)mxGetScalar(pargin[3]);
    }

    mwIndex N = (mwIndex)taylordegree(1.0,tol);
    if (nargin == 5){
        N = (mwIndex)mxGetScalar(pargin[4]);
    }
DEBUGPRINT(("declare outputs\n"));
    pargout[0] = mxCreateDoubleMatrix(n,1,mxREAL);
    pargout[1] = mxCreateDoubleMatrix(1,1,mxREAL);
    pargout[2] = mxCreateDoubleMatrix(1,1,mxREAL);
    
    // decode the sparse matrix
    mwIndex* cp = mxGetJc(A);
    mwIndex* ri = mxGetIr(A);
    double* a = mxGetPr(A);
DEBUGPRINT(("get sparsemat ptrs \n"));
    double* y = mxGetPr(pargout[0]);
    double* nsteps = mxGetPr(pargout[1]);
    double* npushes = mxGetPr(pargout[2]);
    
DEBUGPRINT(("Starting call \n"));
    
    if ( N <= 0 ){
        mexErrMsgIdAndTxt("gexpmq_mex:wrongInputParamterN",
                          "gexpmq_mex needs N > 0");
    }
    if ( tol <= 0 || tol > 1){
        mexErrMsgIdAndTxt("gexpmq_mex:wrongInputParamterTol",
                          "gexpmq_mex needs 0 < tol <= 1");
    }
    if ( c >= n ){
        mexErrMsgIdAndTxt("gexpmq_mex:wrongInputParamterC",
                          "gexpmq_mex needs 1 <= c <= n");
    }
    if ( maxsteps <= 0 ){
        mexErrMsgIdAndTxt("gexpmq_mex:wrongInputParamterN",
                          "gexpmq_mex needs maxsteps > 0");
    }
    
DEBUGPRINT(("calling gexpmq\n"));
    gexpmq(n, cp, ri, a, // sparse matrix
           c, N, tol, maxsteps, // parameters
           y, nsteps, npushes);
    
}
