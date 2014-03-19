
experimentname = 'nexpokit_runtime';

datalist = { 'itdk0304-cc', 'dblp-cc', 'flickr-bidir-cc'};
%, 'ljournal-2008', 'twitter-2010', 'com-friendster'};
num_data = numel(datalist);
num_trials = 5;
num_algs = 6;
heading = '\n \t expmv \t half \t gsqres \t gexpmq \t gexpm_hash \t expm_svec';


maxnnz = 10000;
tol = 1e-4;
t = 1;

datasizes = zeros(num_data,1); % holds n for each graph
seeds = zeros(num_trials,num_data); % holds seeds for each trial for each graph
xtrues = sparse(1,num_trials,num_data);
time_vals = zeros(num_algs,num_trials,num_data);
err_vals = zeros(num_algs,num_trials,num_data);


P=1;

for dataindex = 1:num_data
	clear P;
	dataset = datalist(dataindex);
	load(strcat('/scratch2/dgleich/kyle/colstochdata/', dataset));
	datasizes(dataindex) = size(P,1);
	n = datasizes(dataindex);
	seeds(:,dataindex) = randi(n,num_trials,1);
	if n > size(xtrues,1),
		xtrues(end:n,1)=0;
	end

	fprintf('\n Num of trials = %i \t dataset = %s \t tol = %5.8f \t maxnnz = %i', ...
		num_trials, dataset, tol, maxnnz);
	fprintf(heading);

	for trial=1:num_trials
		fprintf('\n\n trial num = %i',trial);
		ind = seeds(trial);
		ec = eyei(n,ind);

		alg_num = 1;
			tic; [x_true,s,m,mv,mvd] = expmv(t,P,ec,[],'single'); time_vals(alg_num, trial, dataindex) = toc;
			normtrue = norm(x_true,1); %[vtrue strue] = sort(x_true, 'descend');
			xtrues(:,trial,dataindex) = x_true;
		fprintf('\n time = %8.7f', time_vals(trial,alg_num));

		alg_num = alg_num + 1;
			tic; [y,s,m,mv,mvd] = expmv(t,P,ec,[],'half'); time_vals(alg_num, trial, dataindex) = toc;
			err_vals(alg_num, trial, dataindex) = norm(x_true - y,1)/normtrue;
		fprintf('\n time = %8.7f', time_vals(alg_num, trial, dataindex));

		alg_num = alg_num + 1;
			tic; [y npush] = gsqres_mex(P,ind,tol,t); time_vals(alg_num, trial, dataindex) = toc;
			err_vals(alg_num, trial, dataindex) = norm(x_true - y,1)/normtrue;
		fprintf('\t %8.7f', time_vals(alg_num, trial, dataindex));
		
		alg_num = alg_num + 1;
			tic; [y nstep npush] = gexpmq_mex(P,ind,tol,t); time_vals(alg_num, trial, dataindex) = toc;
			err_vals(alg_num, trial, dataindex) = norm(x_true - y,1)/normtrue;
		fprintf('\t %8.7f', time_vals(alg_num, trial, dataindex));

		alg_num = alg_num + 1;
			tic; [y hpush hstep] = gexpm_hash_mex(P,ind,tol,t); time_vals(alg_num, trial, dataindex) = toc;
			err_vals(alg_num, trial, dataindex) = norm(x_true - y,1)/normtrue;
		fprintf('\t %8.7f', time_vals(alg_num, trial, dataindex));

		alg_num = alg_num + 1;
			tic; [y svpush] = expm_svec_mex(P,ind,tol,t,maxnnz); time_vals(alg_num, trial, dataindex) = toc;
			err_vals(alg_num, trial, dataindex) = norm(x_true - y,1)/normtrue;
		fprintf('\t %8.7f', time_vals(alg_num, trial, dataindex));

		fprintf('\n err = '); disp(	err_vals(2:5,trial,dataindex)' );

	end % end trials for that dataset
end % end datasets
	aveerrors = sum(err_vals(:,:,dataindex))/num_trials;
	avetimes = sum(time_vals(:,:,dataindex))/num_trials;
disp(strcat('\n',heading))
fprintf('\n ave error'); disp(aveerrors)
fprintf('\n ave times'); disp(avetimes)	

save(['~/nexpokit/results/' experimentname '.mat'], 'seeds', 'err_vals', 'time_vals', ...
		'datasizes', 'xtrues', '-v7.3');