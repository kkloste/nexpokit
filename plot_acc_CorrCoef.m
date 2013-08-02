%% Load the data

load test_acc_CorrCoef

%%
cc100gexmp = squeeze(record(:,5,:));
boxplot(cc100gexmp','positions',log(squeeze(record(:,1,1))));