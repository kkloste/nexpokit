%%
%
%       ACCURACY: Kendall-tau
%
%   Compare Kendall-tau of top 25 and top 100 nodes vs network size
%
% Measuring performance of NEXPOKIT: rankings vs network size,
% compared with EXPOKIT's expv and mexpv functions
% as well as my implementation of the
% horner's rule method for taylor polynomial approximation of expm

       % how to use RECORD:
       %
       %  record(network_id, attribute_id, trial_number)
       %    network_id = the size of the network (which also tells us
       %                     which network these trials were on)
       %    trial_number = which trial, out of 100, this data point is
       %
       %    attribute_id = determines which attribute you look at,
       %                    runtime, or accuracy, and for which function:
       %        1 - size of network being worked in column network_id
       %        2 - time for 'kmatexp' to compute
       %        3 - gex time
       %        4 - gex CC25
       %        5 - gex CC100
       %        6 - expv time
       %        7 - expv CC25
       %        8 - expv CC100
       %        9 - mexpv time
       %        10 - mexpv CC25
       %        11 - mexpv CC100


record = zeros(15,11,100);


for network_id=1:15

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

    n = size(P,1)
    max_iter = 30*ceil(sqrt(n));
    N = 11;
    tol = 1e-5;

    record(network_id,1,:) = n;

    trials = 100;


    for trial_number=1:trials
        j = randi(n);
    	ej = zeros(n,1);
    	ej(j) = 1;

       % 'TRUTH' (high order taylor approx via horner's rule)
       % index in record(:,2,:)
        tic;
        e = kmatexp(P,j);
        record(network_id,2,trial_number) = toc;
        %record(network_id,2,trial_number) = TIME_horn
        [~,s] = sort(e,'descend');

       % NEXPOKIT
       % index in record(:,[3,4,5],:)
       %    3 - time
       %    4 - Ktau25
       %    5 - Ktau100
        tic;
        x = gexpm_mex(P,j,N,tol,max_iter);
        record(network_id,3,trial_number) = toc;    % TIME_gex

        [~,st] = sort(x, 'descend'); % Ktau25,Ktau100_gex
        record(network_id,4,trial_number) = corr(s(1:25),st(1:25),'type','Kendall');
        record(network_id,5,trial_number) = corr(s(1:100),st(1:100),'type','Kendall');

       % EXPOKIT - expv
       % index in record(:,[6,7,8],:)
       %    6 - time
       %    7 - Ktau25
       %    8 - Ktau100
        tic;
        x = expv(1,P,ej, 1e-5, 15);
        record(network_id,6,trial_number) = toc;    % TIME_expv

        [~,st] = sort(x, 'descend'); % Ktau25,Ktau100_expv
        record(network_id,7,trial_number) = corr(s(1:25),st(1:25),'type','Kendall');
        record(network_id,8,trial_number) = corr(s(1:100),st(1:100),'type','Kendall');


       % EXPOKIT - mexpv
       % index in record(:,[9,10,11],:)
       %    9 - time
       %    10 - Ktau25
       %    11 - Ktau100
        tic;
        x = mexpv(1,P,ej, 1e-5, 15);
        record(network_id,9,trial_number) = toc;    % TIME_mexpv

        [~,st] = sort(x, 'descend'); % Ktau25,Ktau100_mexpv
        record(network_id,10,trial_number) = corr(s(1:25),st(1:25),'type','Kendall');
        record(network_id,11,trial_number) = corr(s(1:100),st(1:100),'type','Kendall');

    end % done all trials for that network

end % done all networks
%%