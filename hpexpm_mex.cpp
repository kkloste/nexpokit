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

/**
 * @param n - sparse matrix size
 * @param cp - sparse matrix column pointer (CSC)
 * @param ari - sparse matrix row index (CSC)
 * @param a - sparse matrix values (CSC)
 * @param c - column index
 * @param N - number of terms of Taylor polynomial
 * @param maxnnz - the maximum number of nnz from input vec to use (1 <= maxnnz <= n)
 * @param y - the output vector (length n)
 */

void hpexpm(const mwSize n, const mwIndex* cp, const mwIndex* ari, const double* a, 
           const mwIndex c, const double N, const mwIndex maxnnz, double* y)
{
    // Note: N is taken to be of type "double" because
    // it's used in computations below as type "double".
    
    mwSize maxlistsize = (mwSize)ceil(100*sqrt(n)*log(maxnnz)/log(10));
	double valind,minval;
    mwIndex hsize = 0;
    mwIndex ind, Tind, ri, arinzi, listind, Nk, nzi;
    mwSize listsize = 0;
    
    // allocate data 

    std::vector<mwIndex> Tvec(maxnnz,0);
    std::vector<double> Vvec(maxnnz,0);
    mwIndex *T = &Tvec[0];
    double *v = &Vvec[0];
	
    std::vector<mwIndex> listvec(maxlistsize,0);
    mwIndex *list = &listvec[0];
    
// mexPrintf("pass initialize \n");
    // Unroll first iterate of Horner's Rule

    for (nzi=cp[c]; nzi < cp[c+1]; ++nzi) {
		arinzi = ari[nzi];
        y[arinzi] += a[nzi]/N;
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
 * Stop as soon as every entry of vec has been checked, or we've added
 * enough entries to make the 1-norm > tol, whichever occurs first.
 *
 * ---- I've temporarily taken out the 'tol' feature until I can debug
 * ----     - Kyle, 10/10/13
 *
 *
 * STEP 2: Do heap-matvec y = A*v
 *
 *  - A is compressed sparse column (CSC)
 *  - for each entry in the heap, add v[T[i]]*A[:,T[i]] to y
 *  - do this for each Horner's Rule iterate, and then we're done.
 *
 ***/
    
    
    
    // HORNER'S RULE ITERATES
    for (mwIndex k=1; k <= N-1 ; k++){
		Nk = (mwIndex)N-k;
        // y = y./(N-k);
        // y = A*y;
		// y[c] = y[c] + 1;
                    listind = 0; hsize = 0; valind = 0.; minval = 0.; ri = 0;
                    // for each entry of y, check if we put it in the heap
                    while ( listind < listsize ){
						ind = list[listind];
                        valind = y[ind];
                        
                        // if the value is nonzero, we might add it to heap
                        if (valind > 0) {
                            // if the heap is full, we might have to delete
                            if ( hsize >= maxnnz ){
                                minval = y[T[0]];
                                // if new entry is big enough, add it
                                if ( minval < valind ){
                                    // replace T[0] with y[ind], then heap_down
                                    ri = T[0];
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
                    for (ind = 0; ind < hsize ; ind++){
                        Tind = T[ind];
                        v[ind] = y[Tind];
                        y[Tind] = 0;
                    }

                 // At this point, y[] = 0 entirely, and v contains
                 // just the entries from y that survived the heap process,
                 // i.e. the maxnnz largest entries that were in y.
                    
					listsize = 0;
                    for (ind = 0; ind < hsize; ind++) {
                        Tind = T[ind];
                        valind = v[ind]/(Nk);
                            for ( nzi=cp[Tind]; nzi < cp[Tind+1]; ++nzi) {
								arinzi = ari[nzi];
								if (y[arinzi] == 0){
									list[listsize] = arinzi;
									listsize++;
								}
                                y[arinzi] += a[nzi]*valind;
                            }
                        v[ind] = 0;
                    }
        if (y[c] == 0){
			list[listsize] = c;
			listsize++;
		}
        y[c] = y[c] + 1;
		
//		mexPrintf("iteration \n");
    }//terms of Taylor are complete
return;
}

void mexFunction(
  int nargout, mxArray *pargout[],       // these are your outputs
  int nargin, const mxArray *pargin[])   // these are your arguments
{
    // inputs
    // A - sparse n-by-n
    // c - positive integer, the column index of exp(A) to be computed
    // N - positive integer, number of terms of Taylor series used
    // maxnnz - integer scalar max number of entries to use in the matvecs 
    
    const mxArray* A = pargin[0];
    mwIndex c = (mwIndex)mxGetScalar(pargin[1])-1;
    double N = (double)mxGetScalar(pargin[2]);
    mwIndex maxnnz = (mwIndex)mxGetScalar(pargin[3]);
    mwSize n = mxGetM(A);
	
    pargout[0] = mxCreateDoubleMatrix(n,1,mxREAL);
    
    // decode the sparse matrix
    mwIndex* cp = mxGetJc(A);
    mwIndex* ri = mxGetIr(A);
    double* a = mxGetPr(A);

    double* y = mxGetPr(pargout[0]);
    
    // mexPrintf("Starting call \n");
    mxAssert( nargin == 5, "incorrect number of inputs");
    mxAssert( N > 0, "N must be bigger than 0");
	mxAssert( 0 <= c &&  c < n, "c must be 1 <= c <= n");
    mxAssert( 1 <= maxnnz && maxnnz <= n, "we must have 1 <= maxnnz <= n");
    
    hpexpm(n, cp, ri, a, // sparse matrix
          c, N, maxnnz, // parameters
          y); // output product vector   
}
