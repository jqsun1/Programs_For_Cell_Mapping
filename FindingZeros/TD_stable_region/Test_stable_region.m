clear all
close all
clc
% -------------------------------------------------------------------------
% Test of the generalized stablility code for 2D and 3D cases.
% -------------------------------------------------------------------------
% [kp, kd]
wq = cputime;
%
lb = [-8;-5];
ub = [8;1];
N_old = [10; 10];

% stb = @(v, tau)Test_PD_Stability(v, tau);
% stb = @(v, tau)Test_PD_Stability_tracking(v, tau);
% stb = @(v, tau)dof2CTApd(v, tau);
% stb = @(v, tau)dof2CTApid(v, tau);

eps = 1; delta = 4; tau = pi/4; k = 20;
stb = @(x,tau)SD_mathieu(x,delta,eps,tau,k)-1;

iter = 3;
%
% hybrid cell mapping approach for stable region finding
% [S, N, h] = Stable_region(stb, lb, ub, N_old, tau);
[S_new, N_new, h] = Stable_region_boundary(stb, lb, ub, N_old, tau, iter);
%
% rule out artificial cells at the final stage (resembles the dominance check)
cell_rmv = [];
for i = 1:length(S_new)
    z = celltoz(S_new(i), N_new);
    x = ztox(z, h, lb);
    MaxEig = stb(x, tau); % CTA stability test
    if MaxEig > 0
        cell_rmv = [cell_rmv; S_new(i)];
    end
end
S_new = setdiff(S_new, cell_rmv); % final ruling out
%
% use PlotCells for final boundary depict
% % estimated lower and upper bound for MOP
% lb = min(x,[],1)' - h/2;
% ub = max(x,[],1)' + h/2;
% [lb, ub]
%
time = cputime - wq
% save FunctionZeros_TD_PD_simplified.mat
% save FunctionZeros_TD_PD_simplified_N_15x15_iter_4.mats