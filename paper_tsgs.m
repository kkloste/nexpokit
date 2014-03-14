function x = gexpm(P,c,tol)
n = size(P,1); N = 11; 
r = zeros(n,N+1); r(c,1) = 1; x = zeros(n,1); sumr=1;  % the residual and solution
while sumr >= tol % use max iteration too
  [ml,q] = max(r(:)); i=mod(q-1,n)+1; k=ceil(q/n);% use a heap in practice for max
  r(q) = 0; x(i) = x(i)+ml; sumr = sumr-ml;% zero the residual, add to solution
  [nset,~,vals] = find(P(:,i)); ml=ml/k;   % look up the neighbors of node i
  for j=1:numel(nset)                                  % for all neighbors
    if k==N, x(nset(j)) = x(nset(j)) + vals(j)*ml;     % add to solution
    else, r(nset(j),k+1) = r(nset(j),k+1) + vals(j)*ml;% or add to next residual
          sumr = sumr + vals(j)*ml;                 
end, end, end                              % end if, end for, end while
