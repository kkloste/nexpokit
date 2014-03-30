dataset = 'webbase-2001';
load /scratch2/dgleich/kyle/data/webbase-2001
A = Problem.A;
clear Problem;
P = normout(A)';
clear A;

save(['/scratch2/dgleich/kyle/colstochdata/' dataset '.mat'],'P','-v7.3');
exit