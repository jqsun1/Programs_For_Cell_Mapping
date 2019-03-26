function [P, scm] = gcm_rand(dsys, lb, ub, N, samNo, samTra)
% -------------------------------------------------------------------------
% Build GCM with random sampling method, the output is the sparse matrix
% representation of the graph. Note the P(i,j) entry is the transition
% probability from cell i to j. Compatibale SCM is also extracted from GCM,
% in the format of a 1d array. This gcm construction code is for random
% system where each sampling point is integrated with a certain number of
% trajectories.
%
% Output arguments:
%       P:   sparse matrix representation of gcm
%      scm:  array representation of extracted scm
%
% By: Free Xiong; 2014-10-13
% -------------------------------------------------------------------------
h = (ub-lb)./N;
rows = [];
cols = [];
p = [];
% cimg = zeros(samNo*samTra,1);
%
Nc = prod(N); % no partition involved, cell range is from 1 to Nc
scm = zeros(Nc,1);
for cs = 1:Nc
    cimg = zeros(samNo,samTra);
    xc = ztox(celltoz(cs,N),h,lb); % use only one inital condition
    lb_cs = xc-h/2;
    ub_cs = xc+h/2;
    
    % random sampling within each cell
    for i = 1:samNo
        x0 = lb_cs + rand(length(h),1).*(ub_cs-lb_cs);
        % generate several trajectories in each sampling point
        for j = 1:samTra
            ximg = dsys(x0);
            cimg(i,j) = ztocell(xtoz(ximg,h,lb),N);
        end
    end
    
    cimg = cimg(:);
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