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

#include "mex.h"
#include <queue>
#include <vector>
#include <assert.h>
#include <math.h>

/** A replacement for std::queue<int> using a circular buffer array */
class array_queue {
    public:
    size_t max_size;
    std::vector<int> array;
    size_t head, tail;
    size_t cursize;
    array_queue(size_t _max_size)
    : max_size(_max_size), array(_max_size), head(0), tail(0), cursize(0)
    {}
    
    void empty() {
        head = 0;
        tail = 0;
        cursize = 0;
    }
    
    size_t size() {
        return cursize;
    }
    
    void push(int i) {
        assert(size() < max_size);
        array[tail] = i;
        tail ++;
        if (tail == max_size) {
            tail = 0;
        }
        cursize ++;
    }
    
    int front() {
        assert(size() > 0);
        return array[head];
    }
    
    void pop() {
        assert(size() > 0);
        head ++;
        if (head == max_size) {
            head = 0;
        }
        cursize --;
    }
};


/**
 * Computes the degree N for the Taylor polynomial
 * of exp(tP) to have error less than eps*exp(t)
 *
 * ( so exp(-t(I-P)) has error less than eps 
 *  or exp(tP) has relative error less than eps
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
 * @param c - the column index to approximate (0 <= c < n)
 * @param tol - the stopping tolerance (0 < tol < Inf)
 * @param y - the output vector (length n)
 * @param nsteps - the number of output steps (length 1)
 * @param npushes - the number of total pushes (length 1)
 */
void gsqres(const mwSize n, const mwIndex* cp, const mwIndex* ari, const double* a,
            const mwIndex c, const double tol,
            double* y, double *nsteps, double *npushes)
{
    //mexPrintf("Input n=%i N=%i c=%i tol=%i maxsteps=%i\n", n, N, c, tol, maxsteps);
    double t = 1; // if we want to allow for parameter t in future, just make it a user-input parameter
    mwIndex N = taylordegree(t,tol);
    std::vector<double> psivec(N+1,0.);
    psivec[N] = 1;
    for (mwIndex k = 1; k <= N ; k++){
        psivec[N-k] = psivec[N-k+1]*t/(double)(N-k+1) + 1;
    } // psivec[k] = psi_k(t)
    double toln = tol/(double)n;
    for (mwIndex k = 0; k <= N ; k++){
        psivec[k] = toln/psivec[k];
    } // psivec[k] = toln/psivec[k]
    
    mwIndex M = n*N;
    
    // allocate data
    std::vector<double> rvec(M,0.);
    double *r = &rvec[0];
    
    std::queue<mwIndex> Q;
    
    // i is the node index, j is the "step"
    #define rentry(i,j) ((i)+(j)*n)
    
    // mexPrintf("Init...\n");
 
    // set the initial residual, add to the heap, and update
    r[rentry(c,0)] = 1;
    Q.push(rentry(c,0));
    
    mwIndex npush = 0;
    mwIndex iter = 0;
    
    //mexPrintf("Loop...\n");
    
    while (1) {
        
        /* STEP 1: pop top element off of queue
         *  * get indices i,j from T
         *  * add r(i,j) to y(i)
         *  * set r(i,j) to zero (update sumresid)
         * STEP 2: get i^th column from A
         *  * get neighbors of ith node
         *  * (if j == N-1), add the column to y instead of r.
         *  * add as a column to next time-step of r, and update queue
         *  *  (update sumresid)
         * Check for convegence!
        */
        
        // STEP 1: pop top element off of heap
        mwIndex ri = Q.front();
        Q.pop();
        
        // decode incides i,j
        mwIndex i = ri%n;
        mwIndex j = ri/n;
        
        double rij = r[ri];
        
        // update yi, zero out residual entry
        y[i] += rij;
        r[ri] = 0;
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
            // this is the interior case, and so we add the column of A 
            // to the residual at the next time step.
            for (mwIndex nzi=cp[i]; nzi < cp[i+1]; ++nzi) {
                mwIndex v = ari[nzi];
                double ajv = a[nzi];
                mwIndex re = rentry(v,j+1);
                double reold = r[re];
                r[re] += ajv*rijs;
                
                if (r[re] >= psivec[j+1]) {
                    if (reold < psivec[j+1]) {
                        Q.push(re);
                    }
                }
            }
            npush+=cp[i+1]-cp[i];
        }
        iter ++;
        if ( Q.size() == 0) {
            *nsteps = (double)iter;
            break;
        }
    }// end "while 1"
    
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
    const mxArray* A = pargin[0];
    mwSize n = mxGetM(A);
    mwIndex c = (mwIndex)mxGetScalar(pargin[1])-1;
    double tol = mxGetScalar(pargin[2]);

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
    
    mxAssert(tol > 0 && tol <= 1, "tol must be 0 < tol <= 1");
    mxAssert(c >= 0 && c < n, "column c must be 1 <= c <= n");
    
    gsqres(n, cp, ri, a, // sparse matrix
           c, tol,  // parameters
           y, nsteps, npushes);
    
}
    
    
  