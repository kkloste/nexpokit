% nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nojvm -nosplash -singleCompThread -r maxnnz_plot > /scratch2/dgleich/kyle/joblog/accmaxplot.txt &
addpath('../plotting_utils'); % so "set_figure_size.m" is available
experimentname = 'acc_maxnnz';
experiment_directory = '/scratch2/dgleich/kyle/nexpokit/results/';

addpath('~/nexpokit/plotting');

load(strcat(experiment_directory, experimentname, '_to_plot') );

% For each value of topks = [25,100,1000]
% try each pair: {1-norm, kendall, intersect} vs {maxnnz, maxnnz/edgedensity}

[num_metrics, num_data, num_maxnnzs, num_topks] = size(acc_data);
acc_xaxis = zeros(num_maxnnzs,1);
maxnnzs = [100,200,500,1000,2000,5000,10000];
topks = [25,100,1000];

figdims = [3,3];

%	%% PLOT : log10(1-norm / edgedensity) (metricsid =1)
%	clf;
%	hold all
%	metricsid = 1; % use 1-norm
%	topksid = 1; % doesn't matter, since 1-norm is same for all
%	for dataid=1:num_data
%		xaxis = log10(maxnnzs);
%	%    xaxis = log10(maxnnzs./edgedensity(dataid));
%	%    xaxis = maxnnzs;
%	%    xaxis = maxnnzs./edgedensity(dataid);
%		plot(xaxis, squeeze(acc_data(metricsid,dataid,:,topksid)) - log10(edgedensity(dataid)) );
%	end
%
%	ylim([-5.5,-1]);
%	title('Error (1-norm) vs. Nonzeros used');
%	xlabel('log10(maxnnz)');
%	ylabel('log10(error / edge-density)');
%	legend('A','B','C','D','E','F','G','Location','Northeast');
%	legend boxoff;
%	set_figure_size(figdims);
%	print(gcf,strcat('maxnnz_vs_norm_edense','.eps'),'-depsc2');


%% PLOT : log10(1-norm) (metricsid =1)
clf;
hold all
metricsid = 1; % use 1-norm
topksid = 1; % doesn't matter, since 1-norm is same for all
for dataid=1:num_data
    xaxis = log10(maxnnzs);
%    xaxis = log10(maxnnzs./edgedensity(dataid));
%    xaxis = maxnnzs;
%    xaxis = maxnnzs./edgedensity(dataid);
    plot(xaxis, squeeze(acc_data(metricsid,dataid,:,topksid)) );
end

ylim([-5.5,-1]);
title('Error (1-norm) vs. Nonzeros used');
xlabel('log10(maxnnz)');
ylabel('log10(error / edge-density)');
legend('A','B','C','D','E','F','G','Location','Northeast');
legend boxoff;
set_figure_size(figdims);
print(gcf,strcat('maxnnz_vs_norm','.eps'),'-depsc2');


%	%% PLOT : kendall (metricsid = 2)
%	clf;
%	hold all
%	metricsid = 2; % usekendall
%	topksid = 3; % 25,100,1000
%	xlineorig = topks(topksid);
%	for dataid=1:num_data
%		xaxis = log10(maxnnzs); xline = log10(xlineorig);
%	%    xaxis = log10(maxnnzs./edgedensity(dataid)); xline = log10(xlineorig);
%	%    xaxis = maxnnzs; xline = xlineorig;
%	%    xaxis = maxnnzs./edgedensity(dataid); xline = xlineorig;
%		plot(xaxis, squeeze(acc_data(metricsid,dataid,:,topksid)));    
%	end
%	hx = graph2d.constantline(xline, 'LineStyle','--','Color',[.7 .7 .7]);
%%
%	title('Kendall Tau accuracy vs. Nonzeros used');
%	xlabel('log10(maxnnz)');
%	ylabel('Kendall Tau');
%	legend('A','B','C','D','E','F','G','Location','Southeast');
%	legend boxoff;
%	set_figure_size(figdims);
%	print(gcf,strcat('maxnnz_vs_Kendall','.eps'),'-depsc2');



%% PLOT : intersect (metricsid = 3)
clf;
hold all
metricsid = 3; % use intersect
topksid = 3; % 25,100,1000
xlineorig = topks(topksid);
% ylim([0,1]);
for dataid=1:num_data
   xaxis = log10(maxnnzs); xline = log10(xlineorig); 
%    xaxis = log10(maxnnzs./edgedensity(dataid)); xline = log10(xlineorig);
%    xaxis = maxnnzs; xline = xlineorig;
%    xaxis = maxnnzs./edgedensity(dataid); xline = xlineorig;
    plot(xaxis, squeeze(acc_data(metricsid,dataid,:,topksid)));    
end
hx = graph2d.constantline(xline, 'LineStyle','--','Color',[.7 .7 .7]);
changedependvar(hx,'x');

title('Set Precision vs. Nonzeros used');
xlabel('log10(maxnnz)');
ylabel('Set Precision');
legend('A','B','C','D','E','F','G','Location','Southeast');
legend boxoff;
set_figure_size(figdims);
print(gcf,strcat('maxnnz_vs_Precision','.eps'),'-depsc2');


%% PLOT : intersect w edge density (metricsid = 3)
clf;
hold all
metricsid = 3; % use intersect
topksid = 3; % 25,100,1000
xlineorig = topks(topksid);
% ylim([0,1]);
for dataid=1:num_data
%   xaxis = log10(maxnnzs); xline = log10(xlineorig); 
    xaxis = log10(maxnnzs./edgedensity(dataid)); xline = log10(xlineorig);
%    xaxis = maxnnzs; xline = xlineorig;
%    xaxis = maxnnzs./edgedensity(dataid); xline = xlineorig;
    plot(xaxis, squeeze(acc_data(metricsid,dataid,:,topksid)));    
end
hx = graph2d.constantline(xline, 'LineStyle','--','Color',[.7 .7 .7]);
changedependvar(hx,'x');

title('Set Precision vs. Nonzeros used');
xlabel('log10(maxnnz/edgedensity)');
ylabel('Set Precision');
legend('A','B','C','D','E','F','G','Location','Southeast');
legend boxoff;
set_figure_size(figdims);
print(gcf,strcat('maxnnz_vs_Precision_edense','.eps'),'-depsc2');

exit