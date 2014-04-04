
%% fix the data
load scaling_1
results = results(1:17,:)';
maxd = full([gdata.maxdeg]);
maxd2 = maxd.^2;
ns = [gdata.n];
nnzs = [gdata.nnz];
gsize = ns+nnzs;
p10 = prctile(results,10);
p25 = prctile(results,25);
p75 = prctile(results,75);
p50 = prctile(results,50);
p90 = prctile(results,90);
%%
clf; hold on;
%patch([gsize fliplr(gsize)],[p10 fliplr(p90)],[0.8,1,0.8]);
%patch([gsize fliplr(gsize)],[p25 fliplr(p75)],[0.5,1,0.5]);
patch([gsize fliplr(gsize)],[p25 fliplr(p75)],[0.7,1,0.7]);

set(gca,'XScale','log');
set(gca,'YScale','log');
plot(gsize,p50,'.-','MarkerSize',18,'LineWidth',2,'Color',[0,0,0.8]);
xlim([min(gsize),max(gsize)]);
set(gca,'XTick',[1e4,1e6,1e8]);
ylabel('time in seconds');
xlabel('graph size');
set_figure_size([3,3]);
print(gcf,'ff-040-runtime.eps','-depsc2');

%%
load scaling_s
results = results(1:15,:)';
gdata = gdata(1:15);
maxd = full([gdata.maxdeg]);
maxd2 = maxd.^2;
ns = [gdata.n];
nnzs = [gdata.nnz];
gsize = ns+nnzs;
p10 = prctile(results,10);
p25 = prctile(results,25);
p75 = prctile(results,75);
p50 = prctile(results,50);
p90 = prctile(results,90);

%%

clf; hold on;
%patch([gsize fliplr(gsize)],[p10 fliplr(p90)],[0.8,1,0.8]);
patch([gsize fliplr(gsize)],[p25 fliplr(p75)],[0.7,1,0.7]);

set(gca,'XScale','log');
set(gca,'YScale','log');
plot(gsize,p50,'.-','MarkerSize',18,'LineWidth',2,'Color',[0.0,0,0.8]);
xlim([min(gsize),max(gsize)]);
ylabel('time in seconds');
xlabel('graph size');
set(gca,'XTick',[1e4,1e6,1e8]);
set_figure_size([3,3]);
print(gcf,'ff-048-runtime.eps','-depsc2');

%%
clf; hold on;
patch([maxd2 fliplr(maxd2)],[p25 fliplr(p75)],[0.7,1,0.7]);

set(gca,'XScale','log');
set(gca,'YScale','log');
plot(maxd2,p50,'.-','MarkerSize',18,'LineWidth',2,'Color',[0.0,0,0.8]);
xlim([min(maxd2),max(maxd2)]);
ylabel('time in seconds');
xlabel('graph size');
set(gca,'XTick',[1e4,1e6,1e8]);
set_figure_size([3,3]);

ylabel('time in seconds');
xlabel('max-degree squared');
line([1e3,1e9],[1e-4,1e2]);

print(gcf,'ff-048-scaling-d2.eps','-depsc2');
