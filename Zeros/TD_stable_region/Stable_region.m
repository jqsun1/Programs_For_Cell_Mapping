function [S, N, h] = Stable_region(stb, lb, ub, N, tau)
% -------------------------------------------------------------------------
% This code implements the stability domain searching by employing the
% hybrid searching pattern. Very coarse grids are firstly used to roughtly
% identify those stable regions, refinements are then conducted to find
% finer stable regions untile certain convergence criteria is met. This
% approach is expected to save a lot of computation effects than SD method.
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
%
% Output arguments:
%       S:   cell set with stability condition statisfied
%       N:   final partition after subdivision
%       h:   final cell size after subdivision
% -------------------------------------------------------------------------
%
h = (ub - lb)./N;
S = (1:prod(N))'; % initial cell set
%
stable_volume_old = prod(ub-lb);
flag = true;
%
counter = 0;
while flag
    cell_rmv = [];
    for i = 1:length(S)
        z = celltoz(S(i), N);
        x = ztox(z, h, lb);
        MaxEig = stb(x, tau); % CTA stability test
        if MaxEig > 0
            cell_rmv = [cell_rmv; S(i)];
        end
    end
    %
    % new set for stable region
    stable_volume_new = prod(h)*(length(S)-length(cell_rmv));
    counter = counter + 1;
    %
    if (stable_volume_old - stable_volume_new)/stable_volume_old < 0.01...
            || counter > 2
        flag = false;
    else
        div = ones(length(N),1); % 2n+1 subdivision
        h = h./(2*div+1);
        % subdivision of the previous stable region and update new
        % partition for the next iteration
        [S, N] = refine(setdiff(S,cell_rmv), lb, ub, N, div);
    end
    %
    stable_volume_old = stable_volume_new;
end