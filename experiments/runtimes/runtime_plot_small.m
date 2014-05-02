
experimentname = 'runtime_small';
addpath('~/nexpokit');
addpath('~/nexpokit/plotting');
load(strcat('~/nexpokit/results/' , experimentname));
addpath('../plotting_utils'); % so "set_figure_size.m" is available
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

%% Now plot

[num_graphs, ~, num_algs] = size(percdata);
% percdata = zeros(num_graphs,3,num_algs);

% if you want to omit a dataset, omit it from 'newindexing'
newindexing = [1:num_graphs];

clf;
hold all;
hs = [];
%alglist = { 'expmv', 'half', 'gexpmq', 'gexpm', 'expmimv'};
colors = 'bcgrmk';
for id=1:num_algs
%	plotstr = [colors(id) 'o'];
%	scatter(log10(inputsize),percdata(:,2,id),plotstr);
%	scatter(log10(inputsize),percdata(:,1,id),[colors(id) '.']);
%	scatter(log10(inputsize),percdata(:,3,id),[colors(id) '.']);

	% guarantee we plot from smallest to largest by input size
	subset = inputsize(newindexing);
	subperc = percdata(newindexing,:,:);
	[~,perm] = sort(subset);
hs(end+1) = plot(log10(subset(perm)), log10(subperc(perm,2,id)),[colors(id) '.-']);
%	plot(log10(subset(perm)), log10(subperc(perm,1,id)),[colors(id) '.--']);
%	plot(log10(subset(perm)), log10(subperc(perm,3,id)),[colors(id) '.--']);	
	%hs(end+1) = plot(log10(subset(perm)), subperc(perm,2,id),[colors(id) '.-']);
	%hs(end+1) = plot(log10(subset(perm)), subperc(perm,1,id),[colors(id) '.--']);
	%hs(end+1) = plot(log10(subset(perm)), subperc(perm,3,id),[colors(id) '.--']);
end

title('Runtime v. Inputsize for tol = 1e-4');
xlabel('log10(|V|+|E|)');
ylabel('log10(runtime) (s)');
legend(hs,'expmv', 'half', 'gexpmq', 'gexpm', 'expmimv','Location','Southeast');
legend boxoff;
%xlim([4.5,9.75]);
set_figure_size([5,3]);
print(gcf,strcat(experimentname,'.eps'),'-depsc2');

exit