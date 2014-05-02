addpath('../plotting_utils'); % so "set_figure_size.m" is available

P = load_staged_data('ljournal-2008');
%%
y = kmatexp(P,1);
n = size(P,1);
%% 
rng(1);
p = randperm(n);
yp = y(p);
%%
clf;
plot(yp,'LineWidth',1.5);
set_figure_size([3,2.5]);
xlim([1,n]);

%%
xlabel(sprintf('nnz = %i', nnz(y)));
ylabel('magnitude');
box off;
%%
print('localize-ljournal.eps','-depsc2');


%%
clf;
[~,p] = sort(y,'descend');
%loglog(y(p)/sum(y),'LineWidth',1.5);
loglog(1-cumsum(y(p))/sum(y(p)),'k-','LineWidth',1.5);
set(gca,'GridLineStyle','-');
xlim([1,n]);
set(gca,'XTick',[1 10 100 1000 10000 100000 1e6]);
set(gca,'YTick',[1e-14 1e-12 1e-10 1e-8 1e-6 1e-4 1e-2 1]);
grid on;
grid minor
ylabel('1-norm error');
xlabel('largest non-zeros retained');

%%
set_figure_size([3,2.5]);
greygrid(gca,[0.55,0.55,0.55]);
%greygrid;
print('localize-ljournal-order.eps','-depsc2');

%%
tols = logspace(-1,-5);
errs1 = zeros(length(tols),1);
work1 = zeros(length(tols),1);
nnzs1 = steps;
for i=1:length(tols)   
    [xapprox npushes nsteps] = gexpm_mex(P,1,tols(i));
    work1(i) = npushes;
    nnzs1(i) = nnz(xapprox);
    errs1(i) = norm(xapprox/sum(xapprox)-y/sum(y),1);
end
for i=1:length(tols)
    [xapprox,npushes] = gexpmq_mex(P,1,tols(i));
    work2(i) = npushes;
    nnzs2(i) = nnz(xapprox);
    errs2(i) = norm(xapprox/sum(xapprox)-y/sum(y),1);
end
%%
nnzsseq = logspace(0.5,4);
for i=1:numel(nnzsseq)
    [xapprox npushes] = expmimv_mex(P,1,1e-5,1.,nnzsseq(i)); 
    work3(i) = npushes;
    nnzs3(i) = nnz(xapprox);
    errs3(i) = norm(xapprox/sum(xapprox)-y/sum(y),1);
end

    
%%
clf;
[~,p] = sort(y,'descend');
%loglog(y(p)/sum(y),'LineWidth',1.5);
hs = [];
hs(end+1) = loglog(1-cumsum(y(p))/sum(y(p)),'.','MarkerSize',2);
hold on;
hs(end+1) = loglog(nnzs1,errs1,'r-','LineWidth',1.25);
hs(end+1) = loglog(nnzs2,errs2,'b--','LineWidth',1.25);
hs(end+1) = loglog(nnzs3,errs3,'g-.','LineWidth',1.25);
colors = [0,0,0,
    117,112,179
27,158,119 %
217,95,2
]/256;
subseq = [1:5:50 50];
loglog(nnzs1(subseq),errs1(subseq),'.','MarkerSize',8,'Color','k');
loglog(nnzs2(subseq),errs2(subseq),'.','MarkerSize',8,'Color','k');
loglog(nnzs3(subseq),errs3(subseq),'.','MarkerSize',8,'Color','k');
for i=1:numel(hs), set(hs(i), 'Color', colors(i,:)); end
set(gca,'GridLineStyle','-');
grid on;
grid minor
xlim([1,n]);
set(gca,'XTick',[1 10 100 1000 10000 100000 1e6]);
set(gca,'YTick',[1e-8 1e-6 1e-4 1e-2 1]);
ylim([1e-8,1]);
legend(hs(2:4),'gexpm','gexpmq','expmimv');
xlabel('non-zeros');
ylabel('1-norm error');
set_figure_size([5,2.5]);
greygrid(gca,[0.55,0.55,0.55]);

print('localize-ljournal-alg.eps','-depsc2');


