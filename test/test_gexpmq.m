
A = load_graph('netscience-cc');
P = normout(A)';
Z = expm(P);
tol = 1e-4;
for j=1:size(P,1)
    [y npush] = gexpmq_mex(P,j,tol,1.,0); % this should be full accuracy
    assert(norm(y-Z(:,j),1)/exp(1) <= tol, sprintf('failed on column %i of netscience ', j))
end


%%
A = load_graph('pgp-cc');
P = normout(A)';
n = size(P,1);
tol = 1e-4;
t = 1;
debugflag=0;
for j=1:50
    [x_true,s,m,mv,mvd] = expmv(1,P,eyei(n,j),[],'single');
    [y npush] = gexpmq_mex(P,j,tol,1.,0); % this should be full accuracy
    assert(norm(y-x_true,1)/norm(x_true,1) <= tol, sprintf('failed on column %i of pgp-cc ', j))
end

%%
A = load_graph('pgp-cc');
P = normout(A)';
n = size(P,1);
tol = 1e-4;
t = 1;
maxpush = 1;
debugflag=0;
for j=1:50
    [y npush] = gexpmq_mex(P,j,tol,1.,0,1); % this should be only itself
    y(j) = y(j) + 1;

    assert(nnz(spones(y)) == 1, sprintf('failed on column %i of pgp-cc ', j))
end
