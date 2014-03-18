
A = load_graph('netscience-cc');
P = normout(A)';
Z = expm(P);
tol = 1e-4;
for j=1:size(P,1)
    [y npush] = gsqres_mex(P,j,tol/2,1.,0); % this should be full accuracy
    assert(norm(y-Z(:,j),1) < tol, sprintf('failed on column %i of netscience ', j))
end