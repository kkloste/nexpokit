warning('superceded by test_steps_accuracy_order');
rng(1); % reset the random generator

ntrials = 10;
topks = [10,25,100,1000];
nets = [11:15];
nsteps = 10;
nmults = [2 5 10 15 25 50];

clear record recordnn recordeffmv recordts recordtsnn

for network_id=nets;
    name = graphnames(network_id);
    A = load_graph(name);
    P = sparse(normout(A)');
    n = size(P,1);
    %p = randperm(n); % generate a random list of nodes
    [~,p] = sort(sum(A,2),'descend');
    p(1:100) = []; % remove the first 100 nodes of highest degree
    maxsteps = ceil(logspace(2,log10(10*n),nsteps));
    tols = logspace(-2,-8,nsteps);
    
    for t=1:ntrials
        j = p(t); % get the tth entry in the random list, it's a ranodm node!
        
        xtrue = kmatexp(P,j,ceil(10*log(n)/log(2))); % compute the exact solution
        [~,px] = sort(xtrue,'descend');
        xtruenn = xtrue; % remove neighbors and self
        xtruenn(j) = -Inf;
        xtruenn(logical(A(:,j))) = -Inf;
        [~,pxnn] = sort(xtruenn,'descend');
        nleft = sum(isfinite(xtruenn));
        
        for si=1:nsteps
            ns = maxsteps(si);
            tol = tols(si);
            [xapprox asteps npushes] = gexpmq_mex(P,j,11,tol,ns);
            
            [~,pxa] = sort(xapprox,'descend');
            xapproxnn = xapprox;
            xapproxnn(j) = -Inf;
            xapproxnn(logical(A(:,j))) = -Inf;
            [~,pxann] = sort(xapproxnn,'descend');
            
            recordeffmv(network_id,si,t) = npushes/nnz(A);
            
            for ki=1:numel(topks)
                k = min(topks(ki),n);
                record(network_id,si,t,ki) = numel(intersect(px(1:k),pxa(1:k)))/k;
            end
            
            for ki=1:numel(topks)
                k = min(topks(ki),nleft);
                recordnn(network_id,si,t,ki) = numel(intersect(pxnn(1:k),pxann(1:k)))/k;
            end
        end
        
        for nmi = 1:numel(nmults)
            nterms = nmults(nmi);
            
            xapprox = kmatexp(P,j,nterms);
            
            [~,pxa] = sort(xapprox,'descend');
            xapproxnn = xapprox;
            xapproxnn(j) = -Inf;
            xapproxnn(logical(A(:,j))) = -Inf;
            [~,pxann] = sort(xapproxnn,'descend');
            
            for ki=1:numel(topks)
                k = min(topks(ki),n);
                recordts(network_id,nmi,t,ki) = numel(intersect(px(1:k),pxa(1:k)))/k;
            end
            
            for ki=1:numel(topks)
                k = min(topks(ki),nleft);
                recordtsnn(network_id,nmi,t,ki) = numel(intersect(pxnn(1:k),pxann(1:k)))/k;
            end
        end
    end
end
save 'test_steps_accuracy.mat' record recordnn recordeffmv recordts recordtsnn nsteps nmults ntrials topks nets;