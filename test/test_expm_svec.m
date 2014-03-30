
A = load_graph('netscience-cc');
P = normout(A)';
n = size(P,1);
Z = expm(P);
tol = 1e-4;
t = 1;
debugflag=0;
maxnnz = floor(n);
for j=1:size(P,1)
    [y npush] = expm_svec_mex(P,j,tol,t,maxnnz,debugflag); % this should be full accuracy
    assert(norm(y-Z(:,j),1)/norm(Z(:,j),1) <= tol, sprintf('failed on column %i of netscience ', j))
end