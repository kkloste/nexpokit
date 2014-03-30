data = load('test_steps_accuracy');

%%
t=4;
ni=data.nets(3);


clf; hold all;
semilogx(squeeze(data.recordeffmv(ni,:,t)),squeeze(data.recordtaunn(ni,:,t,:)),'.-');
semilogx(data.nmults,squeeze(data.recordtstaunn(ni,:,t,2)),'.-');
set(gca,'XScale','log');
title(graphnames(ni));
xlim([1e-3,1e2]);
ylim([-0.1,1.1]);

