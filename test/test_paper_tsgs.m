function check_kmatexp
%%

A = load_graph('karate');
P = normout(A)';

for i=1:size(A,1)
    z = kmatexp(P,i,10); % should be fully accuracy
    z2 = paper_tsgs(P,i,1e-11); % should be "accuracy"
    norm(z-z2)
    assert(norm(z-z2,1)<1e-10, ...
        sprintf('failed on column %i of karate',i));
end

%%

