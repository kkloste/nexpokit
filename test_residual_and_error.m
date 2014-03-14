function test_residual_and_error

%%
N = 11;
tol = 1e-7;
j = 1;

A = load_graph('netscience-cc');
P = normout(A)';
ytrue = kmatexp(P,j,30);


[y,iter,errs,resids] = gexpm_error(P,j,tol,N,ytrue);

semilogy([errs' resids' exp(1)*resids']);
legend('error','resid','eresid');


function [y,iter,errs,resids] = gexpm_error(A,c,tol,N,ytrue)
% gexpm(A) 	returns y as an approximation of
%		the c^th column of e^(a*(A-l*I))
%      	when A has nonnegative off-diagonal entries.
%		Gauss-Southwell (GS) iterative method for solving
%		large sparse systems is used to approximate
%		a truncated power series approximation of e^(a*(A-l*I))
%
% example: input: an n x n matrix A of appropriate format
% 		gexpm(A,c,tol,N,a,l) returns y, an n x 1 column
%		such that e^(a*(A-l*I))*e_c = y, approximately,
%		with the residual in the GS iteration having
%		1-norm < tol. 'N' is the degree of the truncated
%		power series that is used to approximate e^A

% Column-wise expm via Gauss-Southwell iteration
% on truncated power series approximation to e^X
%
% ***** CURRENTLY CODED SPECIFICALLY FOR A = G*D^-1
% inputing any other kind of matrix A will probably
% not work well, and may take a long time to converge



l = 0;
a = 1;
n = size(A,1);

%
% 			VARIABLES
%
%	n = dimension of matrix
%	N = length of the truncated power series
%	r = GS residual vector:
%		r(:,j) holds the (j+1)th block
%			of the GS residual vector
%	c = index of the column of expm(A) we want
%	y = the actual solution that gets returned
%
%	T = tree for the heap
%	L = lookup table for the heap
%	SIZE = # of nonzeros in the heap
%
%	rMAX 	= value of largest entry of r
%	i,j	 	= indices of largest entry of r
%	r1norm	= the 1-norm of the residual of the
%				GS linear system (happens to be 
%				the sum of the entries of r)



r = sparse(zeros(n,N));  
%r = zeros(n,N);

% The data structure for r that is best
% depends on how much memory is needed:
% nonsparse is faster for lookups
% when n*N is small enough.
% ( when n ~ 27,000 nonsparse was faster,
% but for n ~ 190,000 sparse was faster).

y = sparse(zeros(n,1));
T = zeros(n*N,1);
L = zeros(n*N,1);

SIZE=0;
a1=a-1;
N1 = N-1;
aN = a/N;

% FIRST ITERATION
y(c) = 1;
r(:,1) = a*A(:,c);

indices = find(A(:,c));
nnz = size(indices,1);

% This FOR loop re-heapifies the recently added entries
for z=1:nnz
	ind = indices(z);
    %INLINE: [SIZE,T,L]=hupdate(ind,SIZE,T,L,r);
					%		%%%  	 HEAP  UPDATE		%%%		%
					% This inserts r(ind) at the end of the heap,
					% then re-max-heapifies.
					% 	ind = index in r of most recent
					%			entry to be updated

                    					% If r(ind) wasn't in the heap,
                           				% add it at the very end.
						if L(ind)==0, SIZE = SIZE+1;T(SIZE)=ind; L(ind)=SIZE;	end
						p=L(ind);
						% Move entry up the heap, if necessary.
						while 1
							if p==1, break; end			% entry is at the top of heap
							tp=T(p); p2=floor(p/2); tp2=T(p2);	% look at its parent
							if r(tp2)>=r(tp), break;	% parent is larger or equal
							else						% parent is smaller, so swap
								T(p2)=tp; L(tp)=p2;	T(p)=tp2; L(tp2)=p;	p=p2;
							end
						end
		
						% Now p is the current location in T
						% (after possibly moving up) of the
						% entry we just added. Next we move
						% the entry down the heap, if necessary.
						k=p;
						while 1
							q=2*k; tk=T(k);
							if q>SIZE, break; end		% End of heap, so no need
														% to move down further.
							if q==SIZE, lc=T(q);		% Only one child, so skip
														% choice between children.
							else
								lc=T(q); rc=T(q+1);		% Pick largest child
								if r(rc)>r(lc),			% Right child is larger
									q=q+1; lc=rc;
								end
							end
							if r(tk)>=r(lc), break;		% k is larger than or equal
														% to both children, so end
							else
								T(k)=lc; L(lc)=k;		% swap
								T(q)=tk; L(tk)=q;
								k=q;
							end
						end
end

r1norm = 1;
iter = 1;

while (r1norm > tol)
%get max value in r and its indices
	t1 = T(1);
    j = ceil(t1/n);
    i = t1 - (j-1)*n;
    rMAX=r(i,j);    
	% INLINE: [SIZE,T,L]=hpop(SIZE,T,L,r); 
						%		%%%  	 HEAP  POP		%%%		%
						% This swaps the largest entry of r
						% with the entry at the end of the heap,
						% then shrinks the heap SIZE by 1, and
						% finaly re-heapifies.
    						if SIZE==0, SIZE_error = SIZE % occurs if SIZE < 0
        						return;end;
							L(t1) = 0; 
							T(1) = T(SIZE); L(T(SIZE)) = 1;%T(SIZE)=t1; %this seems unnecessary
							SIZE = SIZE-1;
							k=1;
							while 1
								q=2*k; tk=T(k);
								if q>SIZE, break; end		% End of heap, so no need
															% to move down further.
								if q==SIZE, lc=T(q);		% Only one child, so skip
															% choice between children.
								else
									lc=T(q); rc=T(q+1);		% Pick largest child
									if r(rc)>r(lc),			% Right child is larger
										q=q+1; lc=rc;
									end
								end
								if r(tk)>=r(lc), break;		% k is larger than or equal
															% to both children, so end
								else
									T(k)=lc; L(lc)=k;		% swap
									T(q)=tk; L(tk)=q;
									k=q;
								end
							end

% actual iteration:
		% STEP 1: delete largest entry of residual
	r(i,j) = 0;	% this step is expensive when r is
				% a sparse data structure because
				% zeroing an entry takes time to delete it

		% STEP 2: add it to solution vector
    y(i) = y(i) + rMAX;

		% STEP 3: add that entry's neighborhood
		%		to next chunk of residual

	if j==N1;
		y = y + A(:,i)*(rMAX*aN);
		r1norm = r1norm - rMAX;
    else
        dummy = rMAX/(j+1);
        r(:,j+1) = r(:,j+1) + A(:,i)*(dummy*a); % update residual vector
    	r1norm = r1norm - dummy*(a1+j);	% update 1-norm of residual vector
		indices = find(A(:,i));
		nnz = size(indices,1);
		% this FOR loop adds each new residual entry to the heap
		for z=1:nnz
			ind = indices(z) + j*n;
            % INLINE: [SIZE,T,L]=hupdate(ind,SIZE,T,L,r);
							%		%%%  	 HEAP  UPDATE		%%%		%
							% This inserts r(ind) at the end of the heap,
							% then re-max-heapifies.
							% 	ind = index in r of most recent
							%			entry to be updated

								if L(ind)==0,					% If r(ind) wasn't in the heap,
									SIZE = SIZE+1;				% add it at the very end.
									T(SIZE)=ind;
									L(ind)=SIZE;
								end
								p=L(ind);
								% Move entry up the heap, if necessary.
								while 1
									if p==1, break; end			% entry is at the top of heap
									tp=T(p);
									p2=floor(p/2);				% look at its parent
									tp2=T(p2);
									if r(tp2)>=r(tp), break;	% parent is larger or equal
									else						% parent is smaller, so swap
										T(p2)=tp; L(tp)=p2;
										T(p)=tp2; L(tp2)=p;
										p=p2;
									end
								end
		
								% Now p is the current location in T
								% (after possibly moving up) of the
								% entry we just added. Next we move
								% the entry down the heap, if necessary.
								k=p;
								while 1
									q=2*k; tk=T(k);
									if q>SIZE, break; end		% End of heap, so no need
																% to move down further.
									if q==SIZE, lc=T(q);		% Only one child, so skip
																% choice between children.
									else
										lc=T(q); rc=T(q+1);		% Pick largest child
										if r(rc)>r(lc),			% Right child is larger
											q=q+1; lc=rc;
										end
									end
									if r(tk)>=r(lc), break;		% k is larger than or equal
																% to both children, so end
									else
										T(k)=lc; L(lc)=k;		% swap
										T(q)=tk; L(tk)=q;
										k=q;
									end
								end
		end % end FOR loop
    end
   errs(iter) = sum(ytrue-y);
   resids(iter) = r1norm;
   iter = iter + 1;
end


