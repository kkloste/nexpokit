
A = load_graph('netscience-cc', '../data');
P = normout(A)';
n = size(P,1);
Z = expm(P);
tol = 1e-3;
t = 1;
debugflag=0;
maxnnz = floor(n/2);
for j=1:50%size(P,1)
%    [x_true,s,m,mv,mvd] = expmv(1,P,eyei(n,j),[],'single');
[y npush] = expm_svec_mex(P,j,tol,t,maxnnz,debugflag); % this should be full accuracy
assert(norm(y-Z(:,j),1)/norm(Z(:,j),1) < tol, sprintf('failed on column %i of netscience ', j))
%assert(norm(y-x_true,1)/norm(x_true,1) < tol, sprintf('failed on column %i of netscience ', j))
end