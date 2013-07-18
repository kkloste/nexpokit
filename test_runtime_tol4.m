%%
%
%       RUNTIME
%
% Measuring performance of NEXPOKIT: runtime vs network size,
% compared with EXPOKIT's expv and mexpv functions
% as well as my implementation of the
% horner's rule method for taylor polynomial approximation of expm



% List of test graphs, by size:
%       NAME                SIZE        #       NOTES
%	four-clusters           27
%	karate                  34
%	dolphins                62
%	netscience-cc           379         1
%	tapir                   1024        2
%	email                   1133        
%	rand-hyper-4000         3483        3 *** bad CC for gex, investigate
%	ca-HepTh-cc             8638        4
%	pgp-cc                  10680       5
%	ca-AstroPh-cc           17903       6
%	marvel-comics-cc        19365       7
%	as-22july06             22963       8
%	rand-ff-25000-0.4       25000       9
%	rand-ff-25000-0.49      25000
%	cond-mat-2003-cc        27519       10
%	email-Enron-cc          33696       11
%	cond-mat-2005-fix-cc    36458       12
%	soc-sign-epinions-cc    119130      13
%	itdk0304-cc				190914      14
%	dblp-cc                 226413      15
%
%
% note: I used only the numbered networks in my trials.


ntrials = 5;
tol = 1e-5;
N = 11;
methods = { {@(P,j) gexpm_mex(P,j,N,tol,max(1e4,30*size(P,1))), 'TSGS'}
            {@(P,j) gexpmq_mex(P,j,N,tol,max(1e4,30*size(P,1))), 'TSGSQ'}
            {@(P,j) expv(1, P, eyei(size(P,1),j), tol, 15), 'EXPV'}
            {@(P,j) mexpv(1, P, eyei(size(P,1),j), tol, 15), 'MEXPV'}
            {@(P,j) kmatexp(P,j), 'HORNER'}};
        
graphs = 1:15;        
record = zeros(numel(graphs), 2+2*numel(methods));        

for gi=1:numel(graphs)
    
    name = graphnames(graphs(gi));
    A = load_graph(name);
    P = sparse(normout(A)');

    % Parameters of gexpm_mex
    n = size(P,1);
    
    times = zeros(ntrials,numel(methods));
    sums = 0.; % store a sum to keep track of the time
    
    fprintf('Running %i trials on graph %s ...\n', ntrials, name);
    
    for t=1:ntrials
        j = randi(n);
        
        for mi=1:numel(methods)
            f = methods{mi}{1};
            tic
            x = f(P,j);
            times(t,mi) = toc;
            sums = sums + x(1);
        end
    end
            
    
    record(gi,:) = [n nnz(A) mean(times) std(times)];

end

%%
save 'test_runtime.mat' methods record graphs tol ntrials N