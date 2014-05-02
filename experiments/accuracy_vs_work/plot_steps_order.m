function plot_steps_accuracy

data = load('steps_accuracy');
addpath('../plotting_utils'); % so "set_figure_size.m" is available

make_figure(data, 9, 4); % 9th trial on lj
make_figure(data, 4, 1); % 4th trial on dblp
make_figure(data, 3, 2); % 3th trial on flickr
make_figure(data, 25, 2); % 25th trial on flickr
make_figure(data, 30, 2); % 30th trial on flickr
make_figure(data, 50, 3); % 50th trial on itdk
make_figure(data, 1, 3); % 50th trial on itdk


function make_figure(data, t, ni)

%%

clf; hold all;
h = semilogx(squeeze(data.recordeffmv(ni,:,t)),squeeze(data.recordtaunn(ni,:,t,:)),'.-');
set(gca,'XScale','log');
title(data.nets{ni});
xlim([1e-3,1e1]);
ylim([-0.1,1.1]);

set(gca,'LineWidth',0.8);
set(gca,'XTick',[1e-2,1e-1,1]);

set(h(1),'LineWidth',1,'Color','k','LineStyle','--','MarkerSize',15);
set(h(2),'LineWidth',2,'Color','r','LineStyle','-','MarkerSize',18);
set(h(3),'LineWidth',1,'Color','k','LineStyle',':','Marker','+');
set(h(4),'LineWidth',1,'Color','k','LineStyle','-.','Marker','*');
xlabel('Effective matrix-vector products');
ylabel('Kendall tau');
legh = legend('@10','@25','@100','@1000','Location','Southeast');
set(legh,'XColor',[0.9,0.9,0.9],'YColor',[0.9,0.9,0.9]);
leglineobjs = findobj(legh,'type','line');
legxd = get(leglineobjs(2),'XData');
mid = mean(legxd);
set(leglineobjs(2:2:end),'XData',[mid-0.1,mid+0.1]);
legend boxoff;

legp = get(legh,'Position');
legp(1) = legp(1) - 0.05;
legp(2) = legp(2) + 0.1;
set(legh,'Position',legp);

for i=1:numel(data.storetols)
    line([data.recordstol(ni,i,t),data.recordstol(ni,i,t)],[-5,5],'LineWidth',1);
    h = text(data.recordstol(ni,i,t),-0.05,sprintf('tol=10^{%i}',round(log10(data.storetols(i)))),...
        'Rotation',90,'VerticalAlignment','bottom');
end

set_figure_size([3,3]);

print(gcf,sprintf('steps-order-%s-%i.eps',data.nets{ni},t),'-depsc2');

