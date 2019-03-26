function [S, N, h] = Stable_region_boundary(stb, lb, ub, N, tau, iter)
% -------------------------------------------------------------------------
% This code implements the stability domain searching by employing the
% hybrid searching pattern. Very coarse grids are firstly used to roughtly
% identify those stable regions, refinements are then conducted to find
% finer stable regions untile certain convergence criteria is met. This
% approach is expected to save a lot of computation effects than SD method.
% Note this code differs from its previous version that it no longer refine
% the entire stable set, but it refines the stable region boundary.
%
% As a test, this code is initially designed for PID stable region finding
% for nonlinear van-de-pol step response control. Later, this code can be
% transplated to a more general fashion.
%
% Input arguments:
%      stb:  stability function with parameter and time delay as input
%    lb, ub: lower and upper bound of a vast region in parameter space
%       N:   coarse initial partition
%       tau: time delay
%      iter: iteration times for subdivision 
%
% Output arguments:
%       S:   cell set with stability condition statisfied (the boundary)
%       N:   final partition after subdivision
%       h:   final cell size after subdivision
% -------------------------------------------------------------------------
%
h = (ub - lb)./N;
S = (1:prod(N))'; % initial cell set
%
stable_strip_old = prod(ub-lb);
flag = true;
%
counter = 0;
while flag
    cell_rmv = [];
    %
    % S is the strip cell set that marks the boundary
    for i = 1:length(S)
        z = celltoz(S(i), N);
        x = ztox(z, h, lb);
        MaxEig = stb(x, tau); % CTA stability test
        if MaxEig > 0
            cell_rmv = [cell_rmv; S(i)];
        end
    end
    %
    % only the outter boundary cells are taken for further refinement
    strip = [];
    stable_set = setdiff(S,cell_rmv); % coarse strip set in stable region
    for j = 1:length(stable_set)
        cs = stable_set(j);
        neigh = adjcells(cs, N);
        %
        for k = 1:length(neigh)
            z = celltoz(neigh(k), N);
            x = ztox(z, h, lb);
            MaxEig = stb(x, tau);
            if MaxEig >= 0
                strip = [strip; cs];
            end
        end
    end
    S = strip;
    %
    % store neighbours of S as strip set to compensate the holes
    S_expand = S;
    for i = 1: length(S)
        cs = S(i);
        neigh = adjcells(cs, N);
        for j = 1: length(neigh)
            z = celltoz(neigh(j), N);
            x = ztox(z, h, lb);
            MaxEig = stb(x, tau);
            %
            if ~ismember(neigh(j),S_expand) && MaxEig >= 0
                S_expand = [S_expand; neigh(j)];
            end
        end
    end
    S_expand = unique(S_expand);
    S = S_expand;
    %
    % update refinement iteration
    stable_strip_new = prod(h)*length(S);
    counter = counter + 1;
    %
    if (stable_strip_old - stable_strip_new)/stable_strip_old < 0.01...
            || counter > iter
        flag = false;
    else
        div = ones(length(N),1); % 2n+1 subdivision
        h = h./(2*div+1);
        % subdivision of the previous stable region and update new
        % partition for the next iteration
        [S, N] = refine(S, lb, ub, N, div);
    end
    %
    stable_strip_old = stable_strip_new;
end