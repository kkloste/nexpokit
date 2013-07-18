data = load('test_steps_accuracy');

%%
t=10;
ni=data.nets(5);


clf; hold all;
semilogx(squeeze(data.recordeffmv(ni,:,t)),squeeze(data.recordnn(ni,:,t,:)),'.-');
semilogx(data.nmults,squeeze(data.recordtsnn(ni,:,t,2)),'.-');
set(gca,'XScale','log');
title(graphnames(ni));
xlim([1e-3,1e2]);
ylim([-0.1,1.1]);

%%
t=3;
ni=data.nets(4);


clf; hold all;
semilogx(squeeze(data.recordeffmv(ni,:,t)),squeeze(data.recordnn(ni,:,t,:)),'.-');
semilogx(data.nmults,squeeze(data.recordtsnn(ni,:,t,2)),'.-');
set(gca,'XScale','log');
title(graphnames(ni));
xlim([1e-3,1e2]);
ylim([-0.1,1.1]);

%%
t=2;
ni=data.nets(2);


clf; hold all;
semilogx(squeeze(data.recordeffmv(ni,:,t)),squeeze(data.recordnn(ni,:,t,:)),'.-');
semilogx(data.nmults,squeeze(data.recordtsnn(ni,:,t,2)),'.-');
set(gca,'XScale','log');
title(graphnames(ni));
xlim([1e-3,1e2]);
ylim([-0.1,1.1]);