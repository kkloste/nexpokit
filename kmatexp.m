%% Kyle Kloster
% 
%
% COLUMN-WISE MATRIX EXPONENTIAL
% 


function [y] = kmatexp(A,i,M)
% kmexp(A)  returns x as the ith column of e^A
%               when A is square.
% example: input: an nxn matrix A
% kmexp(A) returns x, an nx1 column such that e^A*e_i = x.
format long g;
n = size(A,1);
if nargin<3
    M = ceil(3*log(n)/log(2));
end

y = zeros(n,1);
y(i) = 1/(M+1); 
y = A*y;
y(i) = y(i,1) + 1/M;
y = A*y;

for k=2:M
	y = y./(M+1-k);
    y = A*y;
    y(i) = y(i) + 1;
end


    
