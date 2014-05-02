experimentname = 'degruntime';
experiment_directory = '/scratch2/dgleich/kyle/nexpokit/results/';

addpath('../plotting_utils'); % so "set_figure_size.m" is available
addpath('~/nexpokit/plotting');
load(strcat(experiment_directory , experimentname , '_to_plot') );

% errors ( num_algs, num_degrees, num_graphs )
% times ( num_algs, num_degrees, num_graphs )
% graphsizes ( num_graphs, 1 )
% alldegrees ( num_degrees, num_graphs )

[num_algs, num_degrees, num_graphs] = size(errors);

datalist = { 'itdk0304-cc', 'dblp-cc', 'flickr-scc', 'ljournal-2008', 'webbase-2001', 'twitter_rv-scc', 'com-friendster'};

for graphid = 1:num_graphs
   clf;
   hold all;
   hs = [];

   dataname = char(datalist(graphid));

   colors = 'bcgrmk';
   for id=1:num_algs

	   % guarantee we plot from smallest to largest by input size
	   subset = alldegrees(:,graphid);
	   [~,perm] = sort(subset,'ascend');
	   hs(end+1) = plot(log10(subset(perm)), log10(squeeze(times(id,perm,graphid))),[colors(id) '.-']);
   end

   title(strcat(dataname,', tol=10^-4');
   xlabel('log10(degree)');
   ylabel('log10(runtime) (s)');
   legend(hs,'expmv', 'half', 'gexpmq', 'gexpm', 'expmimv','Location','Southeast');
   legend boxoff;
   %xlim([4.5,9.75]);
   set_figure_size([5,3]);
   print(gcf,strcat('degruntimes_', dataname ,'.eps'),'-depsc2');

end % end GRAPHID for loop
	
%exit