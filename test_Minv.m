n = 15;
M = full(sparse(2:n,1:(n-1),-1./(1:n-1)',n,n)) + eye(n);
inv(M)
c = sum(inv(M))

%%
t1 = c(1);
t = t1;
for i=1:6
    t = i*(t-1) 
    sum(1
end

%%
tN = (n)/(n-1)