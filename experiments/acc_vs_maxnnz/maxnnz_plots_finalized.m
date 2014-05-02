% nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nojvm -nosplash -singleCompThread -r maxnnz_plot > /scratch2/dgleich/kyle/joblog/accmaxplot.txt &
addpath('../plotting_utils'); % so "set_figure_size.m" is available
experimentname = 'acc_maxnnz';
experiment_directory = '../results/';

load(strcat(experiment_directory, experimentname, '_to_plot') );

datasets = {'itdk0304', 'dblp', 'flickr', 'ljournal', ...
                'webbase','twitter','friendster'};

% For each value of topks = [25,100,1000]
% try each pair: {1-norm, kendall, intersect} vs {maxnnz, maxnnz/edgedensity}

[num_metrics, num_data, num_maxnnzs, num_topks] = size(acc_data);
acc_xaxis = zeros(num_maxnnzs,1);
maxnnzs = [100,200,500,1000,2000,5000,10000];
topks = [25,100,1000];

figdims = [3,3];

%%
clf;
hold all
metricsid = 1; % use 1-norm
topksid = 1; % doesn't matter, since 1-norm is same for all
for dataid=1:num_data
    xaxis = maxnnzs;
    yaxis = squeeze(10.^acc_data(metricsid,dataid,:,topksid));
    h = loglog(xaxis, yaxis,'.-','MarkerSize',12);
    did = mod(dataid,4)+2;
    text(xaxis(did),yaxis(did),datasets{dataid},...
        'HorizontalAlignment','left','VerticalAlignment','bottom',...
        'Color',get(h,'Color'));
end
set(gca,'XScale','log');
set(gca,'YScale','log');

xlabel('maxnnz');
ylabel('1-norm error');

legend boxoff;
set_figure_size(figdims);

print(gcf,strcat('maxnnz_vs_norm','.eps'),'-depsc2');

%%
clf;
hold all
metricsid = 3; % use top-k precision
topksid = 3; % 25, 100, 1000
for dataid=1:num_data
    xaxis = maxnnzs./edgedensity(dataid);
    yaxis = squeeze(acc_data(metricsid,dataid,:,topksid));
    h = semilogx(xaxis, yaxis,'.-','MarkerSize',12);
    did = mod(dataid,1)+1;
    text(xaxis(did),yaxis(did),datasets{dataid},...
        'HorizontalAlignment','center','VerticalAlignment','top',...
        'Color',get(h,'Color'));
end
set(gca,'XScale','log');
xlim([1,1e3]);

xlabel('maxnnz/edgedensity');
ylabel('top 1000 precision');

legend boxoff;
set_figure_size(figdims);
print(gcf,strcat('maxnnz_vs_Precision_edense-top-1000','.eps'),'-depsc2');

%%
!cp *.eps ~/Dropbox/publications/nexpokit-journal/nexpokit-JIM/figures/
