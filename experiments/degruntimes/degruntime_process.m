% combine the data from the smaller graphs with twitter and friendster
% /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nojvm -nosplash -r runtime_process > /scratch2/dgleich/kyle/joblog/runtime_proc.txt

experimentname = 'degruntime';
experiment_directory = '/scratch2/dgleich/kyle/nexpokit/results/';

addpath('~/nexpokit');
%load(strcat('~/nexpokit/results/' , experimentname));
load(strcat(experiment_directory, experimentname));

% num_graphs
% num_trials 
% num_algs 

% datasizes = zeros(num_graphs,1); % holds n + nnz(P)/2 for each graph
% seeds = zeros(num_trials,num_graphs); % holds seeds for each trial for each graph
% xtrues = sparse(1,num_trials*num_graphs);
% time_vals = zeros(num_algs,num_trials,num_graphs);
% err_vals = zeros(num_algs,num_trials,num_graphs);

num_graphs = size(datasizes,1);

errors = err_vals;
times = time_vals;
graphsizes = datasizes;

numrecords = num_graphs;

  % now load webbase
  	load(strcat(experiment_directory , experimentname, '_webbase'));
  	numrecords = numrecords + 1;
  	newtrialsnum = size(err_vals,2);
  	errors(:,1:newtrialsnum,numrecords) = err_vals(:,:);
  	times(:,1:newtrialsnum,numrecords) = time_vals(:,:);
  	graphsizes(numrecords,1) = datasizes(:);
  	
  % now load twitter
  	load(strcat(experiment_directory , experimentname, '_twitter'));
  	numrecords = numrecords + 1;
  	newtrialsnum = size(err_vals,2);  	
  	errors(:,1:newtrialsnum,numrecords) = err_vals(:,:);
  	times(:,1:newtrialsnum,numrecords) = time_vals(:,:);
  	graphsizes(numrecords,1) = datasizes(:);
  
  % now load friendster
  	load(strcat(experiment_directory , experimentname, '_friendster'));
  	numrecords = numrecords + 1;
  	newtrialsnum = size(err_vals,2);  	
  	errors(:,1:newtrialsnum,numrecords) = err_vals(:,:);
  	times(:,1:newtrialsnum,numrecords) = time_vals(:,:);
  	graphsizes(numrecords,1) = datasizes(:);  	
	
	
% errors ( num_algs, num_trials, num_data )
% times ( num_algs, num_trials, num_data )
% graphsizes ( num_data, 1 )

num_algs = size(errors, 1);
num_trials = size(errors, 2);
num_graphs = size(graphsizes,1);

percdata = zeros(num_graphs,3,num_algs);
inputsize = graphsizes;

for functionid = 1:num_algs
	for graphid=1:num_graphs
		datax = times(functionid,:,graphid)';
		percdata(graphid,:,functionid) = prctile(datax,[25 50 75],1);
		clear datax;
	end
end

% HAVE ALL DATA for this experiment -- save in plotting/

save(strcat(experiment_directory, experimentname, '_to_plot', '.mat'), 'percdata', 'graphsizes', 'inputsize','-v7.3');
	
exit