function [P, scm] = gcm(dsys, lb, ub, N, samNo)
% -------------------------------------------------------------------------
% Build GCM with random sampling method, the output is the sparse matrix
% representation of the graph. Note the P(i,j) entry is the transition
% probability from cell i to j. Compatibale SCM is also extracted from GCM,
% in the format of a 1d array.
%
% Output arguments:
%       P:   sparse matrix representation of gcm
%      scm:  array representation of extracted scm
%
% By: Free Xiong; 2014-08-15
% -------------------------------------------------------------------------
h = (ub-lb)./N;
n = length(h);
rows = [];
cols = [];
p = [];
cimg = zeros(samNo,1);
%
Nc = prod(N); % no partition involved, cell range is from 1 to Nc
scm = zeros(Nc,1);
for cs = 1:Nc
    x = ztox(celltoz(cs,N),h,lb);
    lb_cs = x-h/2;
    ub_cs = x+h/2;
    
    % random sampling within each cell
    for i = 1:samNo
        x0 = lb_cs + rand(n,1).*(ub_cs-lb_cs);
        ximg = dsys(x0);
        cimg(i) = ztocell(xtoz(ximg,h,lb),N);
    end
    
    % calculate transition probability, sink cells are skiped
    cimg(cimg<1 | cimg>Nc) = [];
    if isempty(cimg)
        continue
    end
    count = length(cimg);
    imgs = unique(cimg);
    nEdges = histc(cimg,imgs);
    for i = 1:length(nEdges)
        rows = [rows; cs];
        cols = [cols; imgs(i)];
        p = [p; nEdges(i)/count];
    end
    scm(cs) = imgs(1); % use the 1st image cell to build scm
end
P = sparse(rows,cols,p,Nc,Nc);
% 
% % transform scm from 1d array to sparse matrix
% rows = 1:Nc;
% rows(scm==0) = [];
% scm(scm==0) = [];
% C = sparse(rows,scm,true,Nc,Nc);