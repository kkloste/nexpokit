dataset = 'flickr-bidir-cc';
A = load_graph(dataset, '/scratch2/dgleich/kyle/data');

P = colnormout(A);
clear A;

save(['/scratch2/dgleich/kyle/colstochdata/' dataset '.mat'],'P','-v7.3');
exit