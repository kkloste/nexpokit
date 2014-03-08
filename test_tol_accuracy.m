
tols = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7];
ntrials = 100;
topks = [25,100];
nets = 5:6;

clear record recordnn;

for network_id=nets;

    if network_id==1,   name = 'netscience-cc';          end
    if network_id==2,   name = 'tapir';                  end
    if network_id==3,   name = 'rand-hyper-4000';        end
    if network_id==4,   name = 'ca-HepTh-cc';            end
    if network_id==5,   name = 'pgp-cc';                 end
    if network_id==6,   name = 'ca-AstroPh-cc';          end
    if network_id==7,   name = 'marvel-comics-cc';       end
    if network_id==8,   name = 'as-22july06';            end
    if network_id==9,   name = 'rand-ff-25000-0.4';      end
    if network_id==10,   name = 'cond-mat-2003-cc';      end
    if network_id==11,   name = 'email-Enron-cc';        end
    if network_id==12,   name = 'cond-mat-2005-fix-cc';  end
    if network_id==13,   name = 'soc-sign-epinions-cc';  end
    if network_id==14,   name = 'itdk0304-cc';           end
    if network_id==15,   name = 'dblp-cc';               end

    
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
            xapprox = gexpm_mex(P,j,11,tol,10*n);
            [~,pxa] = sort(xapprox,'descend');
            xapproxnn = xapprox;
            xapproxnn(j) = -Inf;
            xapproxnn(logical(A(:,j))) = -Inf;
            [~,pxann] = sort(xapproxnn,'descend');
            
            for ki=1:numel(topks)
                k = min(topks(ki),n);
                record(network_id,ti,t,ki) = numel(intersect(px(1:k),pxa(1:k)))/k;
            end
            
            for ki=1:numel(topks)
                k = min(topks(ki),nleft);
                recordnn(network_id,ti,t,ki) = numel(intersect(pxnn(1:k),pxann(1:k)))/k;
            end
        end
    end
end
save 'test_tol_accuracy.mat' record recordnn tols ntrials topks nets;