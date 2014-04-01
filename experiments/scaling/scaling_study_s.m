%% study 2
prob = 0.48;
seq = [1,2,5];
minsize = 1000;
maxsize = 2e9/20; % should give about 2B non-zeros
%maxsize = 2e7/20; % should give about 2B non-zeros
scale = minsize;
sizes = [];
while scale <= maxsize
    sizes = [sizes scale*seq];
    scale = scale*10;
end
sizes(sizes > maxsize) = [];
ntrials = 25;
tol=1e-3;
%%
results = zeros(numel(sizes),ntrials);
gdata = struct;
for si = 1:numel(sizes);
    n = sizes(si); 
    A = graph_model('forest-fire',n,'prob',prob,'initial',10);
    gdata(si).n = n;
    gdata(si).nnz = nnz(A);
    gdata(si).maxdeg = max(sum(A,2));
    startdeg = zeros(ntrials,1);
    for ti=1:ntrials
        p = randi([1,n]);
        startdeg(ti) = sum(A(:,p));
        t0=tic;
        [y npush] = gsqres_mex(A,p,tol,1.,0); 
        dt=toc(t0);
        results(si,ti) = dt;
    end
    gdata(si).startdegs = startdeg;
end
