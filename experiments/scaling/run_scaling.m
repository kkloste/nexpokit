

%% Let's look at scaling up to 1B edges.

% With prob 0.49, we get roughly 20 nnz per row. 
% 
% With prob 0.4, we get roughtly 4 nnz per row.
%
% So that's 2B non-zeros. 
% 
% We generate nnz from forest fire models at roughtly 1.1M edges per
% second, so 2B non-zeros takes 30M or so.

scaling_study_1


%% study2
prob = 0.4;
seq = [1,2,5];
minsize = 1000;
maxsize = 2e9/20; % should give about 2B non-zeros
scale = minsize;
sizes = [];
while scale <= maxsize
    sizes = [sizes scale*seq];
    scale = scale*10;
end
sizes(sizes > maxsize) = [];