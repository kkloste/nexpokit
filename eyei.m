function E = eyei(n,I)
% EYEI Compute the ith canonical basis vector 
% 
%   ei = eyei(n,i) returns the vector of length n where all elements are
%   zero except for the ith entry, which is one. This is the ith column of
%   the n-by-n identity matrix so eyei(n,i) = eye(n)*eyei(n,i);
%
%   E = eyei(n,I) returns a matrix of canonical vectors from a vector of
%   indices i. if I is a row-vector, then E is an n-by-length(i) matrix. if
%   I is a colum-vector, then E is a length(i)-by-n matrix.

% David F. Gleich
% Purdue University, 2013

% :2013-07-17: Initial version

if numel(I) == 1
    E = zeros(n,1);
    E(I) = 1;
else
    % NOT SURE if this is the correct dimensions yet
    if min(size(I)) ~= 1
        error('eyei:notVector', ...
            'I must be a vector argument, not a matrix of size %i-by-%i', ...
            size(I,1), size(I,2));
    elseif ndims(I) > 2
        error('eyei:notVector', ...
            'I must be a vector argument, not a %i-way array', ndims(I));
    end
    if size(I,1) > size(I,2) % then it's a column vector
        E = eyei(n,I')';
    else % then it's a row vector
        E = zeros(n,numel(I));
        for i=1:numel(I)
            E(I(i),i) = 1;
        end
    end
end