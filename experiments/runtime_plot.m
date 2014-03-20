experimentname = 'runtime';

addpath('~/nexpokit/plotting');
load(strcat('~/nexpokit/results/' , experimentname , '_to_plot') );

[num_graphs, ~, num_algs] = size(percdata);
% percdata = zeros(num_graphs,3,num_algs);

% if you want to omit a dataset, omit it from 'newindexing'
newindexing = [1:num_graphs];

clf;
hold all;
hs = [];
%alglist = { 'expmv', 'half', 'gsqres', 'gexpmq', 'gexpm', 'expmsvec'};
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

title('Runtimes for tol = 1e-4');
xlabel('log10(|V|+|E|)');
ylabel('Runtime (s)');
legend(hs,'expmv', 'half', 'gsqres', 'gexpmq', 'gexpm_hash', 'expm_svec','Location','Southeast');
legend boxoff;
%xlim([4.5,9.75]);
set_figure_size([5,3]);
print(gcf,strcat('runtimes','.eps'),'-depsc2');

	
%exit