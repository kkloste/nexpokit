
A = load_graph('pgp-cc');
P = normout(A)';
n = size(P,1);
%Z = expm(P);
tol = 1e-4;
t = 1;
debugflag=0;
for j=1:50
    [x_true,s,m,mv,mvd] = expmv(1,P,eyei(n,j),[],'single');
    [y npush nstep] = gexpm_mex(P,j,tol,t,debugflag); % this should be full accuracy
    assert(norm(y-x_true,1)/norm(x_true,1) <= tol, sprintf('failed on column %i of netscience ', j))
end