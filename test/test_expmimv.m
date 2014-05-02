
A = load_graph('netscience-cc');
P = normout(A)';
n = size(P,1);
Z = expm(P);
tol = 1e-4;
t = 1;
debugflag=0;
maxnnz = floor(n);
for j=1:size(P,1)
    [y npush] = expmimv_mex(P,j,tol,t,maxnnz,debugflag); 
    assert(norm(y-Z(:,j),1)/norm(Z(:,j),1) <= tol, sprintf('failed on column %i of netscience ', j))
end

%%
A = load_graph('pgp-cc');
P = normout(A)';
n = size(P,1);
tol = 1e-4;
t = 1;
maxnnz = 1;
debugflag=0;
for j=1:50
    [y npush] = expmimv_mex(P,j,tol,t,1,debugflag); 
    assert(nnz(y) == nnz(P(:,j))+1, sprintf('failed nnz on column %i of netscience ', j))
end