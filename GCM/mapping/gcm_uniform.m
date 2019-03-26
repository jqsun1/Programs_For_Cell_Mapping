function [P, scm] = gcm_uniform(dsys, lb, ub, N, div)
% -------------------------------------------------------------------------
% Build GCM with uniform sampling method, the output is the sparse matrix
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
%
Nc = prod(N); % no partition involved, cell range is from 1 to Nc
scm = zeros(Nc,1);
for cs = 1:Nc
    % uniform sampling within each cell
    [cimg, N_new]= refine(cs,lb,ub,N,div);
    h_new = (ub-lb)./N_new;
    for i = 1:length(cimg)
        x0 = ztox(celltoz(cimg(i),N_new),h_new,lb);
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