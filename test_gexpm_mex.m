addpath('~/dev/matlab-bgl/');
[A,xy] = load_graph('karate');
P = sparse(normout(A)');
y1 = gexpm(P,1,1e-3,10);
ymex = gexpm_mex(P,1,10,1e-3,10000);