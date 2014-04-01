addpath('~/dev/graph_models/');
%%
A = graph_model('forest-fire',1000,'prob',0.49,'initial',10);
nnz(A)
%%
A = graph_model('forest-fire',10000,'prob',0.49,'initial',10);
nnz(A)
%%
A = graph_model('forest-fire',500000,'prob',0.49,'initial',10);
nnz(A)

%% Let's look at scaling up to 1B edges.

% With prob 0.49, we get roughly 20 nnz per row. 
% 
% With prob 0.4, we get roughtly 4 nnz per per.
%
% So that's 2B non-zeros. 
% 
% We generate nnz from forest fire models at roughtly 1.1M edges per
% second, so 2B non-zeros takes 30M or so.
