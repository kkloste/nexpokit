data = load('test_steps_accuracy');

t=1;
ni=data.nets(4);


clf; hold all;
semilogx(squeeze(data.recordeffmv(ni,:,t)),squeeze(data.recordnn(ni,:,t,:)),'.-');
semilogx(data.nmults,squeeze(data.recordtsnn(ni,:,t,:)),'.-');
set(gca,'XScale','log');