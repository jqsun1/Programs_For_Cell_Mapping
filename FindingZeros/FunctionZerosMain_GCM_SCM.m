clear all
close all
clc
% -------------------------------------------------------------------------
% Hybrid idea in dynamical system's attractor finding. A simple point to
% point difference equation is firstly put under test. This file directly
% finds the attraction domain from state space. Several benchmark examples
% are examined as a primary illustration. To construct the local simple
% cell mapping, different schemes are used. Mostly, poincare intersection
% and finite time integration.
%
% GCM/SCM hybrid covering set finding algorithm is newly integrated to test
% the over all performance. Refined SCM is taken for selection and
% subdivsion operations.
% -------------------------------------------------------------------------
wq = cputime;
%
% up to now the most promising combination is (dopt = 1, sopt = 1 or 7)
sopt = 7; % stepsize control (stepsize.m)
dopt = 1; % local iteration control (dynamic_systems.m)
opt = false; % normalization control
%
max_iter = 8; % max subdivision times
Nm = 1; % forward shooting times (for scm imag holding)
N_old = [13; 8];
div = ones(length(N_old),1); % (2n+1)
No = 1;
[lb, ub, dsys, g, f] = ProblemDef(No, N_old, dopt, sopt, opt); % problem definition with discrete dyn
%
%% First stage coarse finding of covering set
S = covering_set_gcm(N_old, lb, ub, div, dsys, 1:prod(N_old), 'stable');

whos S
time_cover = cputime - wq
% S = add_neighbors(S', N_old, lb, ub);

%% Set-oriented approach: subdivision and selection
% % 1, scm imag holding
% % 2, gcm-scm hybrid; 
% % 3, gcm imag holding; 
% % 4, gcm pre-imag holding (best); 
% % 5, scm pre-imag holding (ok);
% selection_mode = 4;
% [S_new, N_new, iter, mapping_new] = selection_subdivision(max_iter, S, N_old, ...
%     Nm, div, dsys, ub, lb, selection_mode, 'stable');

mode = 3; % refined scm (inherit scm-mop)
[S_new, N_new, iter] = selection_subdivision_old(f, max_iter, S, N_old, Nm, div, dsys, g, ub, lb, mode, []);
%%
whos S_new
time = cputime - wq
% [xa, ~] = Stability_check(S_new, N_new, lb, ub, f, dsys);
% Q = block_classification(S_new, N_new);

% save FunctionZeros_high_dim_Pro12_N_10x10x10x10_iter_3.mat
% save FunctionZeros_poly_N_15x15x15x15_iter_3.mat
% save FunctionZeros_poly_N_10x10x10x10_iter_3.mat
% save FunctionZeros_saddle_N_10x10_iter_4.mat
% save FunctionZeros_TD_PD_N_15x15_iter_4.mat
% save FunctionZeros_2dof_TD_PD_N_5x5x5x5_iter_2.mat
% save FunctionZeros_Mechanism_Pro17_7dim_iter_4.mat
% save FunctionZeros_oliver3D_N_10x10x10_iter_2.mat

% save FunctionZeros_2I2O_N_13x8_iter_2.mat
% save FunctionZeros_2I2O_N_13x8_iter_3.mat
% save FunctionZeros_2I2O_N_13x8_iter_4.mat
% save FunctionZeros_2I2O_N_13x8_iter_5.mat
% save FunctionZeros_2I2O_N_13x8_iter_6.mat
% save FunctionZeros_2I2O_N_13x8_iter_7.mat
% save FunctionZeros_2I2O_N_13x8_iter_8.mat
% save FunctionZeros_2I2O_N_13x8_iter_9.mat
% save FunctionZeros_2I2O_N_13x8_iter_10.mat
% save FunctionZeros_Fractional_N_13x13_iter_4.mat

% figure
% subplot(1,2,1)
% PlotCells(S,lb,ub,N_old,'red','EdgeOn');
% axis square
% subplot(1,2,2)
% PlotCells(S_new,lb,ub,N_new,'red','EdgeOff');
% axis square
% 
% print -depsc FunctionZeros_2I2O_N_13x13_RandSCM_7.eps