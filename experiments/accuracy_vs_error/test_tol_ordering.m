
tols = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7];
ntrials = 100;
topks = [25,100];
nets = [5:6 14 15];

clear record recordnn;

for network_id=nets;

    name = graphnames(network_id);
    A = load_graph(name);
    P = sparse(normout(A)');
    n = size(P,1);
    p = randperm(n); % generate a random list of nodes
    
    for t=1:ntrials
        j = p(t); % get the tth entry in the random list, it's a ranodm node!
        
        xtrue = kmatexp(P,j); % compute the exact solution
        [~,px] = sort(xtrue,'descend');
        xtruenn = xtrue; % remove neighbors and self
        xtruenn(j) = -Inf;
        xtruenn(logical(A(:,j))) = -Inf;
        [~,pxnn] = sort(xtruenn,'descend');
        nleft = sum(isfinite(xtruenn));
        
        for ti=1:numel(tols)
            tol = tols(ti);
            xapprox = gexpmq_mex(P,j,11,tol,10*n);
            xapproxnn = xapprox;
            xapproxnn(j) = -Inf;
            xapproxnn(logical(A(:,j))) = -Inf;
            
            for ki=1:numel(topks)
                k = min(topks(ki),n);
                record(network_id,ti,t,ki) = corr(xapprox(px(1:k)),xtrue(px(1:k)),'type','Kendall');
            end
            
            for ki=1:numel(topks)
                k = min(topks(ki),nleft);
                recordnn(network_id,ti,t,ki) = corr(xapproxnn(pxnn(1:k)),xtruenn(pxnn(1:k)),'type','Kendall');
            end
        end
    end
end
save 'test_tol_ordering.mat' record recordnn tols ntrials topks nets;