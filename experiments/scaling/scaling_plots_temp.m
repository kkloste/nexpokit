addpath('../plotting_utils'); % so "set_figure_size.m" is available
%%
load temp_results
results = results(:,1:30);
%%
maxd = full([gdata.maxdeg]);
ns = [gdata.n];
nnzs = [gdata.nnz];
semilogx(maxd,mean(results'))
%%
semilogx(maxd.^2,mean(results'))
%%
semilogx(ns,mean(results')./(maxd.^2.*log(maxd.^2)))
%%
semilogx(ns,mean(results'))
%%
loglog(maxd.*log(maxd).^(3/2),(mean(results'))./maxd)