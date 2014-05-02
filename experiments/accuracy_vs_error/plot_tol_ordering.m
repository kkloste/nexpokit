function plot_tol_ordering
addpath('../plotting_utils'); % so "set_figure_size.m" is available
%%
data = load('test_tol_ordering.mat');

%%

make_figure(data, data.nets(1), 2);
make_figure(data, data.nets(2), 2);

function make_figure(data, ni, ki)

accs = squeeze(data.recordnn(ni,:,:,2));
h = boxplot(accs','labels',num2str([-2:-1:-7]'));
set(h,'LineWidth',1.3);
ylim([-1,1.1]);
xlabel('log10 of residual tolerance');
ylabel(sprintf('Kendall-\\tau of exact top-%i', data.topks(ki)));
box off;

%%
set_figure_size([3,3]);
print(gcf,sprintf('figures/tol-ordering-%s-%i-nn.eps', graphnames(ni), data.topks(ki)), '-depsc2');

