data = load('test_steps_accuracy');

%%
t=17;
ni=data.nets(5);


clf; hold all;
h = semilogx(squeeze(data.recordeffmv(ni,:,t)),squeeze(data.recordnn(ni,:,t,:)),'.-');
%semilogx(data.nmults,squeeze(data.recordtsnn(ni,:,t,2)),'.-');
set(gca,'XScale','log');
title(graphnames(ni));
xlim([1e-3,1e1]);
ylim([-0.1,1.1]);

set(gca,'LineWidth',0.8);
set(gca,'XTick',[1e-2,1e-1,1]);

set(h(1),'LineWidth',1,'Color','k','LineStyle','--');
set(h(2),'LineWidth',2,'Color','r','LineStyle','-','MarkerSize',15);
set(h(3),'LineWidth',1,'Color','k','LineStyle',':','Marker','+');
set(h(4),'LineWidth',1,'Color','k','LineStyle','-.','Marker','*');
xlabel('Effective matrix-vector products');
ylabel('Precision');
legh = legend('@10','@25','@100','@1000','Location','Southwest');
set(legh,'XColor',[0.9,0.9,0.9],'YColor',[0.9,0.9,0.9]);
leglineobjs = findobj(legh,'type','line');
legxd = get(leglineobjs(2),'XData');
mid = mean(legxd);
set(leglineobjs(2:2:end),'XData',[mid-0.1,mid+0.1]);
legend boxoff;

%legp = get(legh,'Position');
%legp(1) = legp(1) - 0.1;
%set(legh,'Position',legp);

for i=1:numel(data.storetols)
    line([data.recordstol(ni,i,t),data.recordstol(ni,i,t)],[-5,5],'LineWidth',1);
    %h = text(data.recordstol(ni,i,t),0.55,sprintf('tol=',round(log10(data.storetols(i)))));
    %h = text(data.recordstol(ni,i,t),0.45,sprintf('10^{%i}',round(log10(data.storetols(i)))));
    h = text(data.recordstol(ni,i,t),0.45,sprintf('tol=10^{%i}',round(log10(data.storetols(i)))),...
        'Rotation',90,'VerticalAlignment','top');
end

set_figure_size([3,3]);

print(gcf,sprintf('figures/steps-accuracy-%s-%i.eps',graphnames(ni),t),'-depsc2');
    

%%
t=17;
ni=data.nets(5);


clf; hold all;
semilogx(squeeze(data.recordeffmv(ni,:,t)),squeeze(data.recordnn(ni,:,t,:)),'.-');
semilogx(data.nmults,squeeze(data.recordtsnn(ni,:,t,2)),'.-');
set(gca,'XScale','log');
title(graphnames(ni));
xlim([1e-3,1e2]);
ylim([-0.1,1.1]);
for i=1:numel(data.storetols)
    line([data.recordstol(ni,i,t),data.recordstol(ni,i,t)],[-5,5]);
end

%%
t=4;
ni=data.nets(3);


clf; hold all;
semilogx(squeeze(data.recordeffmv(ni,:,t)),squeeze(data.recordnn(ni,:,t,:)),'.-');
semilogx(data.nmults,squeeze(data.recordtsnn(ni,:,t,2)),'.-');
set(gca,'XScale','log');
title(graphnames(ni));
xlim([1e-3,1e2]);
ylim([-0.1,1.1]);

%%
t=19;
ni=data.nets(2);


clf; hold all;
semilogx(squeeze(data.recordeffmv(ni,:,t)),squeeze(data.recordnn(ni,:,t,:)),'.-');
semilogx(data.nmults,squeeze(data.recordtsnn(ni,:,t,2)),'.-');
set(gca,'XScale','log');
title(graphnames(ni));
xlim([1e-3,1e2]);
ylim([-0.1,1.1]);