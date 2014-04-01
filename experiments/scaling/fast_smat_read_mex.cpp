/* 
 *  mex -largeArrayDims -O fast_smat_read_mex.cpp
 * 
 * This code was created for the nexpokit experiments
 * in order to test the forest fire scaling code, we needed
 * a slightly more efficient way to read smat files.
 */

#include <stdlib.h>
#include <string.h>

#include <fstream>

#include "mex.h"
#include "matrix.h"


char* load_string_arg(const mxArray* a, int k)
{
    mwSize buflen;
    char *s;
    int status;

    if (mxIsChar(a) != 1) {
        mexErrMsgIdAndTxt("matlab_bgl:invalidMexArgument",
            "argument %i must be a string", k+1);
    }

    /* Input must be a row vector. */
    if (mxGetM(a) != 1) {
        mexErrMsgIdAndTxt("matlab_bgl:invalidMexArgument",
            "argument %i must be a string (and a row vector)", k+1);
    }

    /* Get the length of the input string. */
    buflen = (mxGetM(a) * mxGetN(a)) + 1;

    /* Allocate memory for input and output strings. */
    s = (char*)mxCalloc(buflen, sizeof(char));

    status = mxGetString(a, s, buflen);
    if (status != 0) {
        mexErrMsgIdAndTxt("matlab_bgl:sizeError",
            "insufficient space to copy argument %i to a string", k+1);
    }

    return s;
}


void mexFunction(int nlhs, mxArray *plhs[], 
        int nrhs, const mxArray *prhs[])
{
    if (nrhs < 1) { mexErrMsgTxt("need at least a filename"); }
    char* filename = load_string_arg(prhs[0],0);
    
    std::ifstream f(filename);
    if (f.is_open()) {
        mwIndex m, n, nnz;
        f >> m >> n >> nnz;
        
        // we load the transpose
        plhs[0] = mxCreateSparseLogicalMatrix(n,m,nnz);
        
        //mwLogical *a = mxGetLogicals(plhs[0]);
        memset(mxGetPr(plhs[0]), 1, nnz); // set all the logicals to true
        
        mwIndex* jc = mxGetJc(plhs[0]);
        mwIndex* ir = mxGetIr(plhs[0]);
        
        mwIndex r = 0, c = 0, firstc = 1;
        jc[0] = 0;
        for (mwIndex nzi = 0; nzi < nnz; ++nzi) {
            mwIndex row, col;
            double val;
            f >> row >> col >> val;
            if (row < r) {
                mexErrMsgIdAndTxt("fastsmat:sorting","need sorted input row = %i, read = %i", r, row);
            }
            while (row != r) {
                r ++;
                jc[r] = nzi;
                c = 0; // reset the last column we saw
                firstc = 1; // and the flag
            }
            if (!firstc && col <= c) { // if this isn't the first column
                mexErrMsgIdAndTxt("fastsmat:sorting","need sorted input row = %i, col = %i, read = %i", r, c, col);
            }
            ir[nzi] = col;
            firstc = 0;
            c = col;
        }
        jc[m] = nnz;
    } else {
        mexErrMsgTxt("could not read filename");
    }
}

