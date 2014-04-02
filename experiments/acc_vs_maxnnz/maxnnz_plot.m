experimentname = 'acc_maxnnz';
experiment_directory = '/scratch2/dgleich/kyle/nexpokit/results/';

load(strcat(experiment_directory, experimentname, '_to_plot') );

% For each value of topks = [25,100,1000]
% try each pair: {1-norm, kendall, intersect} vs {maxnnz, maxnnz/edgedensity}

[num_metrics, num_data, num_maxnnzs, num_topks] = size(acc_Data);
acc_xaxis = zeros(num_maxnnzs,1);
maxnnzs = [100,200,500,1000,2000,5000,10000];
topks = [25,100,1000];

%% PLOT : log10(1-norm) (metricsid =1)
hold all
metricsid = 1; % use 1-norm
topksid = 1; % doesn't matter, since 1-norm is same for all
for dataid=1:num_data
    plot(maxnnzs, squeeze(acc_data(metricsid,dataid,:,topksid)));
end

%% PLOT : kendall (metricsid = 2)
hold all
metricsid = 2; % usekendall
topksid = 2; % 25,100,1000
xlineorig = topks(topksid);
ylim([0,1]);
for dataid=1:num_data
%    xaxis = log10(maxnnzs); xline = log10(xlineorig);
    xaxis = log10(maxnnzs./edgedensity(dataid)); xline = log10(xlineorig);
%    xaxis = maxnnzs; xline = xlineorig;
%    xaxis = maxnnzs./edgedensity(dataid); xline = xlineorig;
    plot(xaxis, squeeze(acc_data(metricsid,dataid,:,topksid)));    
end
hx = graph2d.constantline(xline, 'LineStyle','--','Color',[.7 .7 .7]);
changedependvar(hx,'x');


%% PLOT : intersect (metricsid = 3)
hold all
metricsid = 3; % use intersect
topksid = 3; % 25,100,1000
xlineorig = topks(topksid);
ylim([0,1]);
for dataid=1:num_data
%   xaxis = log10(maxnnzs); xline = log10(xlineorig); 
    xaxis = log10(maxnnzs./edgedensity(dataid)); xline = log10(xlineorig);
%    xaxis = maxnnzs; xline = xlineorig;
%    xaxis = maxnnzs./edgedensity(dataid); xline = xlineorig;
    plot(xaxis, squeeze(acc_data(metricsid,dataid,:,topksid)));    
end
hx = graph2d.constantline(xline, 'LineStyle','--','Color',[.7 .7 .7]);
changedependvar(hx,'x');
