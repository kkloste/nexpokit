% /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nojvm -nosplash -r runtime_plot > /scratch2/dgleich/kyle/joblog/runtime_plot.txt

experimentname = 'runtime';
experiment_directory = '../results/';

addpath('../plotting_utils'); % so "set_figure_size.m" is available
%addpath('~/nexpokit/plotting');
load(strcat(experiment_directory, experimentname, '_to_plot') );

[num_graphs, ~, num_algs] = size(percdata);
% percdata = zeros(num_graphs,3,num_algs);

% if you want to omit a dataset, omit it from 'newindexing'
newindexing = [1:num_graphs];

inputsize = [1406134, 1659333, 9884547, 83354774, 1138045345, 1427920369, 3677742636];

clf;
hold all;
hs = [];

colors = 'kmbrg';
for id=1:num_algs
	% guarantee we plot from smallest to largest by input size
	subset = inputsize(newindexing);
	subperc = percdata(newindexing,:,:);
	[~,perm] = sort(subset);
    hs(end+1) = loglog(subset(perm), subperc(perm,2,id),[colors(id) '.'],'LineWidth',0.5,'MarkerSize',12);
end
set(hs(1),'LineStyle','-','LineWidth',0.5);
set(hs(2),'LineStyle','-','LineWidth',0.5);
set(hs(3),'LineStyle','-.','LineWidth',0.5);
set(hs(4),'LineStyle','-.','LineWidth',0.5);
set(hs(5),'LineStyle',':','LineWidth',0.5);
%title('Runtimes for tol = 1e-4');
xlabel('|V|+ nnz(P)');
ylabel('runtime (s)');
legend(hs,'expmv', 'half', 'gexpmq', 'gexpm', 'expmimv','Location','EastOutside');
legend boxoff;
set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'XTick',[1e5 1e6 1e7 1e8 1e9 1e10]);
%xlim([4.5,9.75]);
set_figure_size([5.5,2.5]);
print(gcf,strcat('runtimes-1e-4','.eps'),'-depsc2');
