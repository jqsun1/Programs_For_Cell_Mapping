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
max_iter = 2; % max subdivision times
Nm = 1; % forward shooting times (for scm imag holding)
N_old = [13; 7];
div = ones(length(N_old),1); % (2n+1)
No = 1;
[lb, ub, dsys, g, f] = ProblemDef(No, N_old, dopt, sopt, opt); % problem definition with discrete dyn
%
%% First stage coarse finding of covering set
S = covering_set_gcm(N_old, lb, ub, div, dsys, 1:prod(N_old), []);

whos S
time_cover = cputime - wq

%% Set-oriented approach: subdivision and selection
% 1, scm forward searching
% 2, gcm backward searching

selection_mode = 1;
[S_new, N_new, iter, mapping_new] = selection_subdivision(max_iter, S, N_old,...
    div, dsys, ub, lb, selection_mode, []);

%%
whos S_new
time = cputime - wq

% postprocessing back to continuous searching space
Extract2dZeros
whos point_sol