function [P, scm] = gcm_uncertain(dsys, lb, ub, lb_par, ub_par, N, samNo)
% -------------------------------------------------------------------------
% Build GCM with random sampling method, the output is the sparse matrix
% representation of the graph. Note the P(i,j) entry is the transition
% probability from cell i to j. Compatibale SCM is also extracted from GCM,
% in the format of a 1d array. Note the sampling is conducted in parameter
% space for the uncertainity modelling. Same initial condition at cell
% center is treated with different parameter set.
%
% Output arguments:
%       P:   sparse matrix representation of gcm
%      scm:  array representation of extracted scm
%
% By: Free Xiong; 2014-09-24
% -------------------------------------------------------------------------
h_par = ub_par - lb_par;
h = (ub-lb)./N;
rows = [];
cols = [];
p = [];
cimg = zeros(samNo*samNo,1);
%
Nc = prod(N); % no partition involved, cell range is from 1 to Nc
scm = zeros(Nc,1);
for cs = 1:Nc
    xc = ztox(celltoz(cs,N),h,lb); % use only one inital condition
    lb_cs = xc-h/2;
    ub_cs = xc+h/2;
    
    % random sampling within each cell
    for i = 1:samNo
        par = lb_par + rand(length(h_par),1).*(ub_par-lb_par);
        for j = 1:samNo
            x0 = lb_cs + rand(length(h),1).*(ub_cs-lb_cs);
            ximg = dsys(x0,par);
            cimg(samNo*(i-1)+j) = ztocell(xtoz(ximg,h,lb),N);
        end
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