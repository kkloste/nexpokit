%%
load test_runtime

%%
nmethods = numel(methods);
loglog(record(:,1), record(:,3:(3+nmethods-1)),'.-');
myleg = {};
for i=1:nmethods
    myleg{i} = methods{i}{2};
end
legend(myleg{:},'Location','Northwest');
legend boxoff;

%%
nmethods = numel(methods);
errorbar(record(:,1)*ones(1,nmethods), record(:,3:(3+nmethods-1)), record(:,(2+nmethods+1):end));
set(gca,'XScale','log');
set(gca,'YScale','log');

%%

%%
load test_runtime
record = record([1,5,10,12,15:16],:)
[~,p] = sort(record(:,1)+record(:,2));
nmethods = numel(methods);
h=loglog(record(p,1)+record(p,2), record(p,3:(3+nmethods-1)),'.-');
set(gca,'XScale','log');
set(gca,'YScale','log');
set(gca,'LineWidth',0.8','XColor',0.55*[1,1,1],'YColor',0.5*[1,1,1]);
set(h(1),'Marker','*','LineStyle','--','LineWidth',1);
set(h(2),'LineStyle','-','LineWidth',2,'Color',[0.8,0,0],'MarkerSize',15);
set(h(3),'LineStyle','-.','LineWidth',1,'Marker','o','Color',[0,0.8,0]);
set(h(4),'LineStyle',':','LineWidth',1,'Marker','o');
set(h(5),'LineStyle','--','Marker','+','LineWidth',1);
ylim([1e-5,10]);
xlim([1e3,9e6]);
myleg = {};
for i=1:nmethods
    myleg{i} = methods{i}{2};
end
myleg{5} = 'TAYLOR';
legend(myleg{:},'Location','Southeast');
box off;
grid on;
grid minor;
set(gca,'GridLineStyle','-');
set(gca,'GridLineStyle','-');
xlabel('|E| + |V|','Color','k');
ylabel('Runtime (secs).','Color','k');
set_figure_size([5,2.5]);
print(gcf,'figures/runtime.eps','-depsc2');