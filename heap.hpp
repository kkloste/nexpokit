/**
 * @file heap.hpp
 * A set of functions to manipulate max-heaps as sets of arrays.
 */


/** Move an element up the heap until it hits the correct place.
 * @param j the index to move
 * @param n the size of the heap
 * @param T the heap items
 * @param L the location of the heap items
 * @param d the values of the heap items
 */
mwIndex heap_up(mwIndex j, mwSize n, mwIndex* T, mwIndex *L, double *d) {
    while (1) {
        if (j==0) { break; } /* the element is at the top */
        mwIndex heapj = T[j];
        mwIndex j2 = (j-1)/2;
        mwIndex heapj2 = T[j2];
        if (d[heapj2] > d[heapj]) {
            break; /* the parent is smaller, so stop */
        } else {
            /* the parent is larger, so swap */
            T[j2] = heapj; L[heapj] = j2;
            T[j] = heapj2; L[heapj2] = j;
            j = j2;
        }
    }
    return j;
}

/** Move an element down the heap until it hits the correct place.
 * @param j the index to move
 * @param n the size of the heap
 * @param T the heap items
 * @param L the location of the heap items
 * @param d the values of the heap items
 */
mwIndex heap_down(mwIndex k, mwSize n, mwIndex* T, mwIndex *L, double *d) {
    mwIndex heapk = T[k];
    while (1) {
        mwIndex i=2*(k+1)-1;
        if (i>=n) { break; } /* end of heap */
        if (i<n-1) { 
            /* pick between children (unneeded if i==heap_size-1) */
            mwIndex left=T[i];
            mwIndex right=T[i+1];
            if (d[right] > d[left]) { 
                i=i+1; /* pick the smaller right child */
            }
        }
        if (d[heapk] > d[T[i]]) {
            /* k is smaller than both children, so end */
            break;
        } else {
            T[k] = T[i]; L[T[i]]=k;
            T[i] = heapk; L[heapk] = i;
            k=i;
        }
    }
    return k;
}

void heap_print(mwSize n, mwIndex* T, mwIndex *L, double *d) {
    mwIndex i;
    for (i=0; i<n; i++) {
        mexPrintf("%6i ", T[i]);
    }
    mexPrintf("\n");
    for (i=0; i<n; i++) {
        mexPrintf("%6i ", L[T[i]]);
    }
    mexPrintf("\n");
    for (i=0; i<n; i++) {
        mexPrintf("%6.4f ", d[T[i]]);
    }
    mexPrintf("\n");
        
}
