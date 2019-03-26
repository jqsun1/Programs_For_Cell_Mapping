function V_new = backward_traj(V_old, N, lb, ub, dsys, tf)
% -------------------------------------------------------------------------
% Calculate the backward influential region of V_old under a given inversed
% dynamical system. The backward integration is discreted in cell space.
% All cells that cover the backward trajectory are picked for each initial
% condition cell in V_old.
%
% Input arguements:
%      V_old:    Original initial condition cell set
%   N, lb, ub:   Cell space partition info
%      dsys:     Dynamical system
%       tf:      Integration time over inversed system
%
% Output arguement:
%      V_new:    New cell set generated based on V_old
% -------------------------------------------------------------------------
V_new = cell(length(V_old),1);
%
h = (ub - lb)./N;
for i = 1:length(V_old)
    cs = V_old(i);
    z = celltoz(cs, N);
    xf = ztox(z, h, lb);
    %
    % inversed integration
    [~, x_inv] = ode45(dsys, [tf, 0], xf);
    %
    bkw = [];
    for j = 1:size(x_inv,1)
        z = xtoz(x_inv(j,:)', h, lb);
        cs = ztocell(z, N);
        bkw = [bkw; cs]; % column vector here for cell2mat usage
    end
    bkw = unique(bkw, 'stable');
    V_new{i,1} = bkw;
end
%
% eliminate reptition over V_new cell array
V_new = cell2mat(V_new);
V_new = unique(V_new);
V_new(ismember(V_new, V_old)) = [];