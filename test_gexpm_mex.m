[A,xy] = load_graph('karate');
P = normout(A)';
y1 = gexpm(P,1,1e-3,10);
ymex = gexpm_mex(P,1,1e-3,10,10000);