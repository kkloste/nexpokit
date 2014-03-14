load /scratch2/dgleich/kyle/colstochdata/flickr-bidir-cc; dataset = 'flickr';
% load /scratch2/dgleich/kyle/colstochdata/ljournal-2008; dataset = 'ljournal';
% load /scratch2/dgleich/kyle/colstochdata/twitter-2010; dataset = 'twitter';
% load /scratch2/dgleich/kyle/colstochdata/com-friendster; dataset = 'friendster';


n = size(P,1);
maxnnz = 10000;
tol = 1e-5
t = 1;
num_trials = 5;
num_algs = 5;
seeds = randi(n,num_trials,1);
time_vals = zeros(num_trials,num_algs);
err_vals = zeros(num_trials, num_algs);


fprintf('\n Num of trials = %i \t dataset = %s',num_trials, dataset);
fprintf('\n ave error = expmv \t gsqres \t gexpmq \t hpexpm \t expm_svec \n');
for trial=1:num_trials
fprintf('\n trial num = %i\n',trial);
    ind = seeds(trial);
    ec = eyei(n,ind);

alg_num = 1;
%fprintf('\n alg_num = %i \t trial = %i',alg_num, trial);
    tic; [x_true,s,m,mv,mvd] = expmv(t,P,ec,[],'half'); time_vals(trial,alg_num) = toc;
    normtrue = norm(x_true,1); %[vtrue strue] = sort(x_true, 'descend');
fprintf('\n time = %8.7f', time_vals(trial,alg_num));
        
alg_num = 2;
%fprintf('\n alg_num = %i \t trial = %i',alg_num, trial);
%    tic; [y npush] = gsqres_mex(P,ind,tol,t); time_vals(trial,alg_num) = toc;
%    err_vals(trial,alg_num) = norm(x_true - y,1)/normtrue;
fprintf('\t %8.7f', time_vals(trial,alg_num));
        
alg_num = 3;
%fprintf('\n alg_num = %i \t trial = %i',alg_num, trial);
    tic; [y nstep npush] = gexpmq_mex(P,ind,tol); time_vals(trial,alg_num) = toc;
    err_vals(trial,alg_num) = norm(x_true - y,1)/normtrue;
fprintf('\t %8.7f', time_vals(trial,alg_num));

alg_num = 4;
%fprintf('\n alg_num = %i \t trial = %i',alg_num, trial);
%    tic; [y hpush hstep] = gexpm_hash_mex(P,ind,tol,t); time_vals(trial,alg_num) = toc;
%    err_vals(trial,alg_num) = norm(x_true - y,1)/normtrue;
    tic; y = hpexpm_mex(P,ind,6,maxnnz); time_vals(trial,alg_num) = toc;
    err_vals(trial,alg_num) = norm(x_true - y,1)/normtrue;
fprintf('\t %8.7f', time_vals(trial,alg_num));

alg_num = 5;
%fprintf('\n alg_num = %i \t trial = %i',alg_num, trial);
    tic; [y svpush] = expm_svec_mex(P,ind,tol,t,maxnnz); time_vals(trial,alg_num) = toc;
    err_vals(trial,alg_num) = norm(x_true - y,1)/normtrue;
fprintf('\t %8.7f', time_vals(trial,alg_num));

fprintf('\n err = %8i \t %8.7f \t %8.7f \t %f \t %f ', 0, err_vals(trial,2),err_vals(trial,3),err_vals(trial,4),err_vals(trial,5));

end
%fprintf('\n ave error = expmv \t gsqres \t gexpmq \t gexpm_hash \t expm_svec \n');
fprintf('\n ave error = expmv \t gsqres \t gexpmq \t hpexpm \t expm_svec \n');
fprintf('\n ave error = %8i \t %8.7f \t %8.7f \t %f \t %f ',0, sum(err_vals(:,2))/num_trials, sum(err_vals(:,3))/num_trials, sum(err_vals(:,4))/num_trials, sum(err_vals(:,5))/num_trials);
fprintf('\n ave time = %8f \t %8f \t %8f \t %f \t %8.7f \n',sum(time_vals(:,1))/num_trials, sum(time_vals(:,2))/num_trials, sum(time_vals(:,3))/num_trials, sum(time_vals(:,4))/num_trials, sum(time_vals(:,5))/num_trials );