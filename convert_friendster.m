dataset = 'com-friendster';
% A = load_graph(dataset, '/scratch2/dgleich/kyle/data');
load /scratch/dgleich/friendster/friendster
P = colnormout(A);
clear A;

save(['/scratch2/dgleich/kyle/colstochdata/' dataset '.mat'],'P','-v7.3');
exit