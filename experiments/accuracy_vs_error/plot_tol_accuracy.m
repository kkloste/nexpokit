function plot_tol_accuracy
addpath('../plotting_utils'); % so "set_figure_size.m" is available
%%
data = load('tol_accuracy.mat');

%%

make_figure(data, 1, 1);
make_figure(data, 2, 1);
make_figure(data, 3, 1);
make_figure(data, 4, 1);

make_figure(data, 1, 2);
make_figure(data, 2, 2);
make_figure(data, 3, 2);
make_figure(data, 4, 2);


%make_figure(data, data.nets(2), 2);

function make_figure(data, ni, ki)

accs = squeeze(data.recordnn(ni,:,:,ki));
h = boxplot(accs','labels',num2str([-2:-1:-8]'));
set(h,'LineWidth',1.3);
ylim([-0.1,1.1]);
xlabel('log10 of error tolerance');
ylabel(sprintf('Precision at %i', data.topks(ki)));
box off;

title(data.nets(ni));

%%
set_figure_size([3,3]);
print(gcf,sprintf('tol-vs_accuracy-%s-%i-nn.eps', data.nets{ni}, data.topks(ki)), '-depsc2');

