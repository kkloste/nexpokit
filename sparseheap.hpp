/**
 * @file sparseheap.hpp
 * A set of functions to manipulate a max-heap as a set of 3 hashtables
 */

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


struct sparsemaxheap{
    typedef unsigned int uint;
    typedef tr1ns::unordered_map< uint, double > map_type;
    typedef tr1ns::unordered_map< uint, uint > T_map;
    typedef tr1ns::unordered_map< int, uint > L_map;
    map_type map;
    T_map T;
    L_map L;
    unsigned int hsize;
    
    /** Remove max entry in heap, replace with last entry
     *  then heap down. Sets val = maxval, returns index of maxval
     */
    uint extractmax(double &val){
        // what if hsize = 0 ?
        /*        if (hsize==0){
         return -1;
         }*/
        uint oldind = T[0];
        val = map[oldind];
        map[oldind] = 0;
        T[0] = T[hsize-1]; // move the last entry to the root
        L[T[0]] = 0;        //update L to reflect this
        //            map.erase(T[hsize-1]); // erase that entry from the heap
        //            T.erase(hsize-1);
        //            L.erase(oldind);
        L[oldind] = -1;
        hsize--;
        heap_down(0);   //then heap down to re-heapify.
        return oldind;
    }
    
    void heap_down(uint k){
        uint heapk = T[k];
        while (1) {
            uint i=2*(k+1)-1;
            if (i>=hsize) { break; } /* end of heap */
            if (i<hsize-1) {
                /* pick between children (unneeded if i==heap_size-1) */
                uint left=T[i];
                uint right=T[i+1];
                if (map[right] > map[left]) {
                    i=i+1; /* pick the larger child */
                }
            }
            if (map[heapk] > map[T[i]]) {
                /* k is larger than both children, so end */
                break;
            } else {
                T[k] = T[i]; L[T[i]]=k;
                T[i] = heapk; L[heapk] = i;
                k=i;
            }
        }
    }
    
    /** When altering an entry of map, this checks
     *  whether that entry is already in map, or if a new entry
     *  need be made in the heap before adding 'value' to map[index].
     */
    void update(uint index, double value){
        if ( L.count(index) == 0 || L[index] == -1){
            // if 'index' was never in L, or if it was but has been removed: insert new entry
            map[index] = value;
            T[hsize] = index;
            L[index] = hsize;
            hsize++;
            heap_up(hsize-1);
        }
        else{ // update old entry
            map[index] += value;
            heap_up(L[index]);
        }
    }
    
    uint heap_up(uint j){
        while (1) {
            if (j==0) { break; } /* the element is at the top */
            uint heapj = T[j];
            uint j2 = (j-1)/2;  /* otherwise, compare with its parent */
            uint heapj2 = T[j2];
            if (map[heapj2] > map[heapj]) {
                break; /* the parent is larger, so stop */
            } else {
                /* the parent is smaller, so swap */
                T[j2] = heapj; L[heapj] = j2;
                T[j] = heapj2; L[heapj2] = j;
                j = j2;
            }
        }
        return j;
    }
};