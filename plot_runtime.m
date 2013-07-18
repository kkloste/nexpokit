%%
load test_runtime

%%
nmethods = numel(methods);
loglog(record(:,1), record(:,3:(3+nmethods-1)),'.-');

%%
nmethods = numel(methods);
errorbar(record(:,1)*ones(1,nmethods), record(:,3:(3+nmethods-1)), record(:,(2+nmethods+1):end));
set(gca,'XScale','log');
set(gca,'YScale','log');