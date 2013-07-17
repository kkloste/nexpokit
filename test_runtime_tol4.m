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


record = zeros(15,14);


for m=1:15

    if m==1,   name = 'netscience-cc';          end
    if m==2,   name = 'tapir';                  end
    if m==3,   name = 'rand-hyper-4000';        end
    if m==4,   name = 'ca-HepTh-cc';            end
    if m==5,   name = 'pgp-cc';                 end
    if m==6,   name = 'ca-AstroPh-cc';          end
    if m==7,   name = 'marvel-comics-cc';       end
    if m==8,   name = 'as-22july06';            end
    if m==9,   name = 'rand-ff-25000-0.4';      end
    if m==10,   name = 'cond-mat-2003-cc';      end
    if m==11,   name = 'email-Enron-cc';        end
    if m==12,   name = 'cond-mat-2005-fix-cc';  end
    if m==13,   name = 'soc-sign-epinions-cc';  end
    if m==14,   name = 'itdk0304-cc';           end
    if m==15,   name = 'dblp-cc';               end

    
    A = load_graph(name);
    P = sparse(normout(A)');

    % Parameters of gexpm_mex
    n = size(P,1)
    max_iter = 30*ceil(sqrt(n));    
    N = 3*ceil(log(n)/log(2));
    tol = 1e-5;


    TIME_horn = 0;
    TIME_gex = 0;
    TIME_expv = 0;
    TIME_mexpv = 0;
    CC25_gex = 0;
    CC50_gex = 0;
    CC25_expv = 0;
    CC50_expv = 0;
    CC25_mexpv = 0;
    CC50_mexpv = 0;

    trials = 10;


    for i=1:trials
        j = randi(n);
    	ej = zeros(n,1);
    	ej(j) = 1;

        tic;
        e = kmatexp(P,j);
        TIME_horn = TIME_horn + toc;
        [p,s] = sort(e,'descend');


% N is the degree of the Taylor polynomial approximation
        tic;
        x = gexpm_mex(P,j,N,tol,max_iter);
        TIME_gex = TIME_gex + toc;

        [~,st] = sort(x, 'descend');
        R = corrcoef(s(1:25),st(1:25));
        CC25_gex = CC25_gex + R(1,2);
        R = corrcoef(s(1:50),st(1:50));
        CC50_gex = CC50_gex + R(1,2);

        tic;
        x = expv(1,P,ej, 1e-5, 15);
        TIME_expv = TIME_expv + toc;

        [~,st] = sort(x, 'descend');
        R = corrcoef(s(1:25),st(1:25));
        CC25_expv = CC25_expv + R(1,2);
        R = corrcoef(s(1:50),st(1:50));
        CC50_expv = CC50_expv + R(1,2);

        tic;
        x = mexpv(1,P,ej, 1e-5, 15);
        TIME_mexpv = TIME_mexpv + toc;

        [~,st] = sort(x, 'descend');
        R = corrcoef(s(1:25),st(1:25));
        CC25_mexpv = CC25_mexpv + R(1,2);
        R = corrcoef(s(1:50),st(1:50));
        CC50_mexpv = CC50_mexpv + R(1,2);
    end

    TIME_horn = TIME_horn/trials;

    TIME_gex = TIME_gex/trials;
    CC25_gex = CC25_gex/trials;
    CC50_gex = CC50_gex/trials;

    TIME_expv = TIME_expv/trials;
    CC25_expv = CC25_expv/trials;
    CC50_expv = CC50_expv/trials;

    TIME_mexpv = TIME_mexpv/trials;
    CC25_mexpv = CC25_mexpv/trials;
    CC50_mexpv = CC50_mexpv/trials;

    % record (:,1) = size of matrix
    % record (:,2) = TIME_horn
    % record (:,3:5) = gex
    % record (:,6:8) = expv
    % record (:,9:11) = mexpv
    % record (m,12) = TIME_horn/TIME_gex
    % record (m,13) = TIME_expv/TIME_gex
    % record (m,14) = TIME_mexpv/TIME_gex

    record(m,1) = n;
    record(m,2) = TIME_horn;
    record(m,3) = TIME_gex;
    record(m,4) = CC25_gex;
    record(m,5) = CC50_gex;
    record(m,6) = TIME_expv;
    record(m,7) = CC25_expv;
    record(m,8) = CC50_expv;
    record(m,9) = TIME_mexpv;
    record(m,10) = CC25_mexpv;
    record(m,11) = CC50_mexpv;
    record (m,12) = TIME_horn/TIME_gex;
    record (m,13) = TIME_expv/TIME_gex;
    record (m,14) = TIME_mexpv/TIME_gex;

end

%% Runtime vs Network Size:

hold all
plot( record(:,1), record(:,5) );
plot( record(:,1), record(:,8) );
plot( record(:,1), record(:,11) );


