clear all
close all
clc
% -------------------------------------------------------------------------
% Use graph theory to find GCM PG cells
% -------------------------------------------------------------------------
%
% up to now the most promising combination is (dopt = 1, sopt = 1 or 7)
sopt = 7; % stepsize control (stepsize.m)
dopt = 1; % local iteration control (dynamic_system.m)
opt = false; % normalization control

mapping_mode = 2; % 1, scm; 2, selected scm; 3, gcm
selection_mode = 3; % 1, set intersection; 2, sequence selection; 3, refined scm

max_iter = 5; % max subdivision times
Nm = 20; % forward shooting times

N_old = [17; 17];
% N_old = [10; 10; 10; 10];
div = 2*ones(length(N_old),1); % (2n+1)

[lb, ub, dsys, g, f] = ProblemDef(7, N_old, dopt, sopt, opt); % problem definition with discrete dyn
GCM = gcmDataBase(N_old, lb, ub, div, dsys);
[G, S, C] = gcm2graph(GCM); % G is a DAG

