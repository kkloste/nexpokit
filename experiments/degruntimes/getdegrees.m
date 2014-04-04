% nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nojvm -nosplash -singleCompThread -r getdegrees > /scratch2/dgleich/kyle/joblog/getdeg.txt &

experimentname = 'degrees';

datalist = { 'itdk0304-cc', 'dblp-cc', 'flickr-scc', 'ljournal-2008'}
alglist = { 'expmv', 'half', 'gexpmq', 'gexpm', 'expmimv'};

addpath('~/nexpokit');


num_data = numel(datalist);
num_trials = 2;
num_algs = numel(alglist);

disp(experimentname);
maxnnz = 10000;
tol = 1e-4;
t = 1;

% LOAD PARAMETERS
%	this contains t, tol, maxnnz, num_trials, whichnodes
run('~/nexpokit/experiments/degruntimes/degruntime_parameters');

		
datasizes = zeros(num_data,1); % holds n for each graph
seeds = zeros(num_trials,num_data); % holds seeds for each trial for each graph
degrees = zeros(num_trials,num_data);

P=1;

for dataindex = 1:num_data
	clear P;
	dataset = char(datalist(dataindex));
	load(strcat('/scratch2/dgleich/kyle/colstochdata/', dataset));
	n = size(P,1);
	
	d = zeros(n,1);
	for i=1:n, d(i) = nnz(P(:,i)); end
	[dummy perm] = sort(d,'descend');	
	seeds(:,dataindex) = perm(whichnodes(1:num_trials)); % the seeds are the largest degree nodes
	degrees(:,dataindex) = dummy(whichnodes(1:num_trials));
	
end % end datasets

save(['/scratch2/dgleich/kyle/nexpokit/results/' experimentname  '.mat'], 'degrees', '-v7.3');
		
exit