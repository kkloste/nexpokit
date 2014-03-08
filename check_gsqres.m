A = load_graph('netscience-cc', '../data');
P = normout(A)';
Z = expm(P);
tol = 1e-4;
for j=1:size(P,1)
[y npush] = gsqres_mex(P,j,tol); % this should be full accuracy
assert(norm(y-Z(:,j)) < tol, sprintf('failed on column %i of netscience ', j))
end