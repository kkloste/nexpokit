dataset = 'flickr-scc';
A = load_graph(dataset, '/scratch2/dgleich/kyle/data');

P = normout(A)';
clear A;

format long g
n = size(P,1)
nnz = nnz(P)
nnzOVERn = nnz/n


save(['/scratch2/dgleich/kyle/colstochdata/' dataset '.mat'],'P','-v7.3');
exit