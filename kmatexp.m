function [y] = kmatexp(A,i,M)
% KMATEXP Compute a column of the matrix expontial via Horner's rule
%
% y = kmatexp(A,i) computes e^A*ei using ceil(3*log(n)/log(2)) terms of
% the Taylor series polynomial for exp(x). This is a good procedure when
% the norm of A is bounded by 1 -- as it is for stochastic matrices -- or
% small. 
%
% y = kmatexp(A,i,k) using k terms of the Taylor series instead
%
% Example:
%   A = load_graph('karate');
%   P = normout(A)';
%   y = kmatexp(P,1);
%   z = expm(P)*eyei(size(P,1),1);
%   norm(y-z)

% History


n = size(A,1);
if nargin<3
    M = ceil(3*log(n)/log(2));
end


y = zeros(n,1);
y(i) = 1/(M+1); 
y = A*y;
y(i) = y(i) + 1;

for k=1:M
    y = y./(M+1-k);
    y = A*y;
    y(i) = y(i) + 1;
end


    
