function plot_tol_accuracy

%%
data = load('test_tol_accuracy.mat');

%%

make_figure(data, data.nets(1), 2);
make_figure(data, data.nets(2), 2);

function make_figure(data, ni, ki)

accs = squeeze(data.recordnn(ni,:,:,2));
h = boxplot(accs','labels',num2str([-2:-1:-7]'));
set(h,'LineWidth',1.3);
ylim([-0.1,1.1]);
xlabel('log10 of residual tolerance');
ylabel(sprintf('Precision at %i', data.topks(ki)));
box off;

%%
set_figure_size([3,3]);
print(gcf,sprintf('figures/tol-accuracy-%s-%i-nn.eps', graphnames(ni), data.topks(ki)), '-depsc2');

