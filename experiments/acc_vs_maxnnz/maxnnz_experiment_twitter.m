experimentname = 'acc_maxnnz_twitter';

datalist = {'twitter_rv-scc'};

addpath('~/nexpokit');


disp(experimentname);
maxnnz = [100,200,500,1000,2000,5000,10000];
topks = [25,100,1000];
tol = 1e-4;
t = 1.0;

num_data = numel(datalist);
num_trials = 50;
num_maxnnzs = numel(maxnnz);
num_ks = numel(topks);

% LOAD PARAMETERS
%	this contains t, tol, maxnnz, num_trials
% run('~/nexpokit/experiments/runtimes/runtime_parameters');
	
edgedensity = zeros(num_data,1); % holds nnz(P)/n for each graph
seeds = zeros(num_trials,num_data); % holds seeds for each trial for each graph
degrees = zeros(num_trials,num_data); % holds out-degree of each seed
xtrues = sparse(1,num_trials*num_data);
time_vals = zeros(num_maxnnzs,num_trials,num_data);
err_vals = zeros(num_maxnnzs,num_trials,num_data);

kendall_vals = zeros(num_data,num_trials,num_maxnnzs,num_ks);
kendall_valsnn = zeros(num_data,num_trials,num_maxnnzs,num_ks);
intersect_vals = zeros(num_data,num_trials,num_maxnnzs,num_ks);
intersect_valsnn = zeros(num_data,num_trials,num_maxnnzs,num_ks);

P=1;

for dataindex = 1:num_data
	clear P;
	dataset = char(datalist(dataindex));
	load(strcat('/scratch2/dgleich/kyle/colstochdata/', dataset));
	n = size(P,1);
	edgedensity(dataindex) = nnz(P)/size(P,1);
	seeds(:,dataindex) = randi(n,num_trials,1);
	if n > size(xtrues,1),
		xtrues(end:n,1)=0;
	end

	fprintf('\n\n Dataset = %s \t tol = %7.6f', dataset, tol);

	for trial=1:num_trials
		fprintf('  %i', trial);
		
		ind = seeds(trial,dataindex);
		ec = eyei(n,ind);
		degrees(trial,dataindex) = nnz(P(:,ind));
			[xtrue,s,m,mv,mvd] = expmv(t,P,ec,[],'single');
			normtrue = norm(xtrue,1);
			[vtrue strue] = sort(xtrue, 'descend');
			xtrues(:,trial + num_trials*(dataindex-1)) = sparse(xtrue);

		% get true ordering/accuracy information
        [~,px] = sort(xtrue,'descend'); % with node and neighbors
        xtruenn = xtrue; % remove neighbors and self
        xtruenn(ind) = -Inf;
        xtruenn(logical(P(:,ind))) = -Inf;
        [~,pxnn] = sort(xtruenn,'descend');
        nleft = sum(isfinite(xtruenn));
			
		for maxnnzval = 1:num_maxnnzs,
			tic; [y svpush] = expmimv_mex(P,ind,tol,t,maxnnz(maxnnzval));
			time_vals(maxnnzval, trial, dataindex) = toc;
			err_vals(maxnnzval, trial, dataindex) = norm(xtrue - y,1)/normtrue;
			
			xapproxnn = y;
			xapprox = y;
			% remove node and its neighbors
            xapproxnn(ind) = -Inf;
            xapproxnn(logical(P(:,ind))) = -Inf;
	        [~,pxa] = sort(xapprox,'descend'); % with node and neighbors
	        [~,pxann] = sort(xapproxnn,'descend'); % without node and neighbors
                    
            % KENDALL
            for ki=1:num_ks
                k = min(topks(ki),n);
                kendall_vals(dataindex,trial,maxnnzval,ki) = corr(xapprox(px(1:k)),xtrue(px(1:k)),'type','Kendall');
            end
            
            for ki=1:num_ks
                k = min(topks(ki),nleft);
                kendall_valsnn(dataindex,trial,maxnnzval,ki) = corr(xapproxnn(pxnn(1:k)),xtruenn(pxnn(1:k)),'type','Kendall');
            end
            
            % INTERSECT
            for ki=1:num_ks
                k = min(topks(ki),n);
                intersect_vals(dataindex,trial,maxnnzval,ki) = numel(intersect(px(1:k),pxa(1:k)))/k;
            end
            
            for ki=1:num_ks
                k = min(topks(ki),nleft);
                intersect_valsnn(dataindex,trial,maxnnzval,ki) = numel(intersect(pxnn(1:k),pxann(1:k)))/k;
            end

			
		end % end different values of maxnnz for that one trial
		
	end % end trials for that dataset
end % end datasets

save(['/scratch2/dgleich/kyle/nexpokit/' experimentname '_expmv_vectors' '.mat'], 'xtrues', 'seeds', '-v7.3');
save(['/scratch2/dgleich/kyle/nexpokit/results/' experimentname '.mat'], 'degrees', 'seeds', 'err_vals', 'time_vals', ...
		'edgedensity', 'err_vals', 'kendall_vals', 'kendall_valsnn', ...
		'intersect_vals', 'intersect_valsnn', '-v7.3');
		
exit