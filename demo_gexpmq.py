"""
demo_gexpmq.py
 
A demonstration of the Gauss-Seidel coordinate relaxation scheme
to estimate a column of the matrix exponential that operates
only on significant entries.  
 
Written by Kyle Kloster and David F. Gleich
"""
 
import collections
import sys
import math
 
# setup the graph
G = {
  1:set([ 2, 3, 5, 6,]),
  2:set([ 1, 4,]),
  3:set([ 1, 6, 7,]),
  4:set([ 2, 5, 7, 8,]),
  5:set([ 1, 4, 6, 8, 9, 10,]),
  6:set([ 1, 3, 5, 7,]),
  7:set([ 3, 4, 6, 9,]),
  8:set([ 4, 5, 9,]),
  9:set([ 5, 7, 8, 20,]),
  10:set([ 5, 11, 12, 14, 15,]),
  11:set([ 10, 12, 13, 14,]),
  12:set([ 10, 11, 13, 14, 15,]),
  13:set([ 11, 12, 15,]),
  14:set([ 10, 11, 12, 25,]),
  15:set([ 10, 12, 13,]),
  16:set([ 17, 19, 20, 21, 22,]),
  17:set([ 16, 18, 19, 20,]),
  18:set([ 17, 20, 21, 22,]),
  19:set([ 16, 17,]),
  20:set([ 9, 16, 17, 18,]),
  21:set([ 16, 18,]),
  22:set([ 16, 18, 23,]),
  23:set([ 22, 24, 25, 26, 27,]),
  24:set([ 23, 25, 26, 27,]),
  25:set([ 14, 23, 24, 26, 27,]),
  26:set([ 23, 24, 25,]),
  27:set([ 23, 24, 25,]),
}
Gvol = 102
eps = 1.e-3
 
## Estimate column c of 
## the matrix exponential vector 
# G is the graph as a dictionary-of-sets, 
# eps is set to stopping tolerance
def compute_psis(N):
    psis = {}
    psis[N] = 1.
    for i in xrange(N-1,-1,-1):
        psis[i] = psis[i+1]/(float(i+1.))+1.
    return psis    
def compute_threshs(eps, N, psis):
    threshs = {}
    threshs[0] = (math.exp(1)*eps/float(N))/psis[0]
    for j in xrange(1, N+1):
        threshs[j] = threshs[j-1]*psis[j-1]/psis[j]
    return threshs
## Setup parameters and constants
N = 11  
c = 1 # the column to compute
psis = compute_psis(N)
threshs = compute_threshs(eps,N,psis)
## Initialize variables
x = {} # Store x, r as dictionaries
r = {} # initialize residual
Q = collections.deque() # initialize queue
sumresid = 0.    
r[(c,0)] = 1.
Q.append(c)
sumresid += psis[0]
## Main loop
for j in xrange(0, N):
  qsize = len(Q)
  pushtol = threshs[j]/float(qsize)
  for qi in xrange(0, qsize):
    i = Q.popleft() # v has r[(v,j)] ...
    rij = r[(i,j)]
    if rij < pushtol:
      continue
    # perform the relax step
    if i not in x: x[i] = 0.
    x[i] += rij
    r[(i,j)] = 0.
    sumresid -= rij*psis[j]
    update = (rij/(float(j)+1.))/len(G[i])
    for u in G[i]:   # for neighbors of i
      next = (u, j+1)
      if j == N-1: 
        if u not in x: x[u] = 0.
        x[u] += update
      else:
        if next not in r: 
            r[next] = 0.
            Q.append(u)
        r[next] += update
        sumresid += update*psis[j+1]
      if sumresid < eps: break
  if len(Q) == 0: break
  if sumresid < eps: break
  
for v in xrange(1,len(G)+1):
    if v in x:
        print "x[%2i] = %.16lf"%(v, x[v])
    else:
        print "x[%2i] = -0."%(v)
print 