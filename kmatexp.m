function [y] = kmatexp(A,i,M)
% kmexp(A)  returns x as the ith column of e^A
%               when A is square.
% example: input: an nxn matrix A
% kmexp(A) returns x, an nx1 column such that e^A*e_i = x.
%
% Kyle Kloster
% COLUMN-WISE MATRIX EXPONENTIAL
% uses a degree M Taylor polynomial to approximate e^A*e_i,
% via a Horner's rule upate.
%

n = size(A,1);
if nargin<3
    M = ceil(3*log(n)/log(2));
end

% initialize
y = zeros(n,1);
y(i,1) = 1;

for k=0:(M-1)
    y = y./(M-k);
    y = A*y;
    y(i,1) = y(i,1) + 1;
end


    
