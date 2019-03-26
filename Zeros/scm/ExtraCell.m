function [S, step] = ExtraCell(dsys, div, S0, N, lb, ub)
% -------------------------------------------------------------------------
% Bring in extra cells using sampling method introduced in GCM. The extra
% cells is taken as the sampling image cells of the original attraction
% cells acquired by SCM. Several steps of this GCM-like extra cell brought
% in might be needed. This approach is suggested by Prof. Sun specially in
% dealing with strange attractors in nonlinear dynamics.
%
% Later, this code is modified as the PG-like searching, where the step
% value is indetermined at the very begining. Main idea of this PG-like
% search is to sample and shoot until the (i+1) step cell set is strictly
% contained in the union of all previous steps' cell sets. This very large
% GCM-like set actually has no information of mapping. A refined SCM scheme
% will be used at the second stage try to enhance the accuracy of this set.
%
% Input arguements:
%      dsys:   Discrete point-to-point mapping
%      div:    Sampling refinement (2n+1)
%      S0:     Original attractor acquired by SCM (might be incomplete)
%      N:      Cell space partition (coarse)
%    lb, ub:   Lower and upper bound of computational region
%
% Output arguments:
%      S:      Union of coarse cells ready to be refined for next stage
%     step:   Sampling times (how many steps for GCM-like mapping)
% -------------------------------------------------------------------------
step = 0;
S_old = S0;
h = (ub - lb)./N;
S_temp = [];
flag = true;
%
while flag
    S_new = [];
    S_temp = [S_temp; S_old];
    %
    for i = 1:length(S_old)
        cell = S_old(i);
        [cs, Ndiv] = refine(cell, lb, ub , N, div);
        hdiv = (ub - lb)./Ndiv;
        %
        for j = 1:length(cs)
            csj = cs(j);
            z = celltoz(csj, Ndiv);
            xa = ztox(z, hdiv, lb); % central points within "cell" (sampled)
            xd = dsys(xa);
            z = xtoz(xd, h, lb);
            img = ztocell(z, N); % image cell of xd in coarse structure
            S_new = [S_new; img];
        end
    end
    %
    S_new = unique(S_new);
    if all(ismember(S_new, S_temp)) == 1 % self-closed set
        break
    else
        S_old = S_new;
        step = step + 1;
    end
end
%
S = unique(S_temp);
Nc = prod(N);
S(~ismember(S,1:Nc)) = []; % in case the expansion includes outter cells