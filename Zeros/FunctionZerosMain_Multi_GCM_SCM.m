clear all
close all
clc
% -------------------------------------------------------------------------
% Hybrid GCM-SCM is used with multiple stages of GCM searching for covering
% set. The generalized two stage searching enables the GCM/SCM searching
% level being selected by users. kFF random sampling technique is used for
% GCM construction with fixed sampling set. This might be more effectively
% for higher dimensional problems.
% -------------------------------------------------------------------------
wq = cputime;
%
% up to now the most promising combination is (dopt = 1, sopt = 1 or 7)
sopt = 7; % stepsize control (stepsize.m)
dopt = 1; % local iteration control (dynamic_systems.m)
opt = false; % normalization control
%
k = 15; % fixed number of sampling points in gcm
max_iter_gcm = 1; % max subdivision times for gcm
max_iter_scm = 14; % max subdivision times for scm (rolling cut at each edge, multiple of the dimension)

N_old = [13; 7]; % initial cell partition can be very coarse
div = ones(length(N_old),1); % (2n+1)
No = 1; % Problem No.
[lb, ub, dsys, g, f] = ProblemDef(No, N_old, dopt, sopt, opt); % problem definition with discrete dyn
%
%% First stage coarse finding of covering set
[S, N_mid] = covering_set_gcm_gradual(N_old, lb, ub, div, dsys, k, max_iter_gcm, 1:prod(N_old));
% [S, N_mid] = covering_set_gcm_gradual_uniform(N_old, lb, ub, div, dsys, max_iter_gcm, 1:prod(N_old));
whos S
N_mid'
time_cover = cputime - wq

%% Set-oriented approach: subdivision and selection
if max_iter_scm == 0
    S_new = S;
    N_new = N_mid;
else
    mode = 3; % refined scm (inherit MOP)
%     mode = 4; % directed search with MOP formulation

    % rolling cutting of each edge
    n = length(N_old); % problem dim
    S1 = S;
    N1 = N_mid;
    sol_vol = prod((ub - lb)./N_mid)*length(S); % solution volume in cell space
    sol_perct = length(S)/prod(N_mid); % solution occupation under current partition
    div = diag(div);
    
    for i = 1:max_iter_scm
        row = mod(i,n);
        if row == 0
            row = n; % refine only one dimension each iteration
        end
        [S_new, N_new] = selection_subdivision_old(f, 1, S1, N1, [], div(row,:)', dsys, g, ub, lb, mode, []);
        %
        % record convergence status at each scm iteration
        sol_vol = [sol_vol; prod((ub-lb)./N_new)*length(S_new)];
        sol_perct = [sol_perct; length(S_new)/prod(N_new)];
        %
        % swap data for next iteration
        S1 = S_new;
        N1 = N_new;
    end
end
%%
whos S_new
N_new'
time = cputime - wq

% postprocessing
Extract2dZeros();
whos point_sol

% save FunctionZeros_HighDimPoly_gcmRand2_scmPermute12_N_8x8x8x8x8x8.mat
% save FunctionZeros_2Dmechanism_gcmRand2_scmPermute14_N_10x10x10x10x10x10x10.mat