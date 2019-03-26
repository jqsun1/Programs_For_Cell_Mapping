clear all
close all
clc
wq = cputime;
% -------------------------------------------------------------------------
% Hybrid idea in dynamical system's attractor finding. A simple point to
% point difference equation is firstly put under test. This file directly
% finds the attraction domain from state space. Several benchmark examples
% are examined as a primary illustration. To construct the local simple
% cell mapping, different schemes are used. Mostly, poincare intersection
% and finite time integration.
% -------------------------------------------------------------------------
%
% up to now the most promising combination is (dopt = 1, sopt = 1 or 7)
sopt = 7; % stepsize control (stepsize.m)
dopt = 1; % local iteration control (dynamic_systems.m)
opt = false; % normalization control

mapping_mode = 2; % 1, scm; 2, selected scm; 3, gcm; 4, max and min probability of s_scm
selection_mode = 3; % 1, set intersection; 2, sequence selection; 3, refined scm;
max_iter = 2; % max subdivision times
Nm = 20; % forward shooting times

N_old = [10; 10; 10; 10];
div = ones(length(N_old),1); % (2n+1)

[lb, ub, dsys, g, f] = ProblemDef(12, N_old, dopt, sopt, opt); % problem definition with discrete dyn
%% First stage coarse finding of P-K cells
S = covering_set(mapping_mode, dsys, g, ub, lb, N_old, div);
% [~, V] = attraction_basin(S, N_old, lb, ub, dsys);
% S = [S;setdiff(1:prod(N_old),[S;V])'];

%% First refine and id block identification
[S1, N1, ~] = selection_subdivision(1, S, N_old, Nm, div, dsys, g, ub, lb, selection_mode, []);
Q1 = block_classification(S1, N1); % blocks in N1 grid system

%% Successive refine and selection based on the refined result
if length(Q1) == 1
    [S_new, N_new, iter] = selection_subdivision_old(max_iter, S1, N1, Nm, div, dsys, g, ub, lb, selection_mode, []);
else
    [S_new, N_new, iter] = selection_subdivision_old(max_iter, S1, N1, Nm, div, dsys, g, ub, lb, selection_mode, Q1);
end

%%
time = cputime - wq;
disp(['Total iteration time for subdivision is:',' ',num2str(iter+1)]);
disp(['CPU running time is:',' ',num2str(time),'s']);
FunctionZerosPlot;
% [x, lambda] = Stability_check(S_new, N_new, lb, ub, f, dsys);
%%
% save FunctionZeros_2I2O.mat
% save FunctionZeros_saddle.mat
% save FunctionZeros_TD_PD.mat
% save FunctionZeros_high_dim.mat
% save FunctionZeros_2I2O_N_13x13.mat
% save FunctionZeros_2I2O_N_13x9.mat