dataset = 'twitter_rv-scc';
addpath('~/nexpokit');

A = load_graph(dataset, '/scratch2/dgleich/kyle/data/');
P = colnormout(A);
clear A;

n = size(P,1)
nnz = nnz(P)
nnzOVERn = nnz/n

save(['/scratch2/dgleich/kyle/colstochdata/' dataset '.mat'],'P','-v7.3');
exit