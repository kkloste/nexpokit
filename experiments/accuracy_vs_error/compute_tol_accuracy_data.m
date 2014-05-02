
rng(1); % reset the random generator

tols = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8];
ntrials = 100;
topks = [25,100];
nets = {'dblp-cc','flickr-scc','itdk0304-cc','ljournal-2008'};

clear record recordnn;

for nid=1:numel(nets);
    net = nets{nid};
    P = load_staged_data(net);
    A = P;
    n = size(P,1);
    p = randperm(n); % generate a random list of nodes
    
    fprintf('starting trials on net %s\n',net);
    
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
            xapprox = gexpmq_mex(A,j,tol,1.,0); 
            [~,pxa] = sort(xapprox,'descend');
            xapproxnn = xapprox;
            xapproxnn(j) = -Inf;
            xapproxnn(logical(A(:,j))) = -Inf;
            [~,pxann] = sort(xapproxnn,'descend');
            
            for ki=1:numel(topks)
                k = min(topks(ki),n);
                record(nid,ti,t,ki) = numel(intersect(px(1:k),pxa(1:k)))/k;
            end
            
            for ki=1:numel(topks)
                k = min(topks(ki),nleft);
                recordnn(nid,ti,t,ki) = numel(intersect(pxnn(1:k),pxann(1:k)))/k;
            end
            fprintf('net = %-15s  trial = %3i   tol = %e\n', net, t, tol);
        end
        save 'tol_accuracy.mat' record recordnn tols ntrials topks nets;
        pause(1);
    end
end
save 'tol_accuracy.mat' record recordnn tols ntrials topks nets;
