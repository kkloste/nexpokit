addpath('../plotting_utils'); % so "set_figure_size.m" is available

A = load_graph('flickr-bidir-cc');
%%
P = normout(A)';
y = kmatexp(P,1);

%% 
p = randperm(size(A,1));
yp = y(p);
%%
plot(yp,'LineWidth',1.5);
set_figure_size([3,2.5]);

%%
print('localize-flickr.eps','-depsc2');

%%
clf;
[~,p] = sort(y,'descend');
%loglog(y(p)/sum(y),'LineWidth',1.5);
loglog(1-cumsum(y(p))/sum(y(p)),'LineWidth',1.5);
set(gca,'GridLineStyle','-');
grid on;
grid minor

%%
set_figure_size([3,2.5]);
greygrid(gca,[0.55,0.55,0.55]);
%greygrid;
print('localize-flickr-order.eps','-depsc2');

%%
tols = logspace(-1,-5);
errs = zeros(length(tols),1);
work = zeros(length(tols),1);
steps = zeros(length(tols),1);
nnzs = steps;
for i=1:length(tols)
    [xapprox,asteps,npushes] = gexpmq_mex(P,1,11,tols(i),10*size(P,1));
    work(i) = npushes;
    steps(i) = asteps;
    nnzs(i) = nnz(xapprox);
    errs(i) = norm(xapprox/sum(xapprox)-y/sum(y),1);
end
    
%%
clf;
[~,p] = sort(y,'descend');
%loglog(y(p)/sum(y),'LineWidth',1.5);
loglog(1-cumsum(y(p))/sum(y(p)),'LineWidth',1.5);
hold on;
loglog(nnzs,errs,'r.--','LineWidth',1.5);
set(gca,'GridLineStyle','-');
grid on;
grid minor

%%
set_figure_size([3,2.5]);
greygrid(gca,[0.55,0.55,0.55]);
print('localize-flickr-alg.eps','-depsc2');

%%
clf;
plot(nnzs,work,'.-')
