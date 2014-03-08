% Kyle Kloster
% 
%
% compute the number of entries that the vectors
% x and y, taken as sets, share.
% 


function [number] = setoverlap(x,y)
% setoverlap(x,y)  returns percentage (of total size of x)
%                      of the total number of elements
%                   that x and y (viewed as sets of elements
%                    rather than vectors of entries) have
%                   in common.

n = max(size(x));
m = max(size(y));
assert(m==n, 'dimensions of x and y must match');

incommon = 0;

for k=1:n
    value = x(k);
    xnum = size(find(x==value),1);
    ynum = size(find(y==value),1);
    if min(size(find(y==value)))>0,
        incommon = incommon+min(xnum,ynum)/xnum;
    end
end
number = incommon/n;



    
