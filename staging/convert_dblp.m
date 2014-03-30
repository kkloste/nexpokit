% nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nojvm -nosplash -r convert_dblp > convdblp.txt &

dataset = 'dblp-cc';
A = load_graph(dataset, '/scratch2/dgleich/kyle/data');

P = colnormout(A);
clear A;

save(['/scratch2/dgleich/kyle/colstochdata/' dataset '.mat'],'P','-v7.3');
exit