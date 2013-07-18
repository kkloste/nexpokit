function check_kmatexp

%% Do a quick unit-test of kmatexp for accuracy

alpha=1;
A = load_graph('karate');
P = normout(A)';
Z = expm(alpha*P);

for j=1:size(P,1)
    y = kmatexp(alpha*P,j,20); % this should be fully accuracy
    assert(norm(y-Z(:,j)) < 10*eps(1), ...
        sprintf('failed on column %i of karate (alpha=%f)', j,alpha))
end
    
%%    
alpha=0.5;
A = load_graph('karate');
P = normout(A)';
Z = expm(alpha*P);

for j=1:size(P,1)
    y = kmatexp(alpha*P,j,20); % this should be fully accuracy
    assert(norm(y-Z(:,j)) < 10*eps(1), ...
        sprintf('failed on column %i of karate (alpha=%f)', j,alpha))
end
   
%%    
alpha=1.25;
A = load_graph('karate');
P = normout(A)';
Z = expm(alpha*P);

for j=1:size(P,1)
    y = kmatexp(alpha*P,j,20); % this should be fully accuracy
    assert(norm(y-Z(:,j)) < 10*eps(1), ...
        sprintf('failed on column %i of karate (alpha=%f)', j,alpha))
end

%%
alpha=1;
A = load_graph('netscience-cc');
P = normout(A)';
Z = expm(alpha*P);

for j=1:size(P,1)
    y = kmatexp(alpha*P,j,20); % this should be fully accuracy
    assert(norm(y-Z(:,j)) < 10*eps(1), ...
        sprintf('failed on column %i of netscience (alpha=%f)', j,alpha))
end