clear all
close all
clc
% -------------------------------------------------------------------------
% Hybrid idea with adaptive cell sizes allowed here.
% -------------------------------------------------------------------------
wq = cputime;
%
N_old = [11; 11];
[lb, ub, dsys] = ProblemDef(3, N_old); % problem definition with discrete dyn
%
%% First stage coarse finding of P-K cells
[~, C, ~] = graph(dsys, ub, lb, N_old);
[Gr, Pe, St] = search(C, N_old);
scm_att = find(Gr~=1 & St==0);

%% Expand the P-K cell set into a complete PG-like set
div = ones(length(N_old),1); % (2n+1)
[S, ~] = ExtraCell(dsys, div, scm_att, N_old, lb, ub); % include extra cells to reveal complete info

%% Set-oriented approach: subdivision and selection (adaptive)
iter = 3;
[Q, P] = Adaptive_attractor(S, N_old, lb, ub, dsys, iter, div);

time = cputime - wq

save DuffingAttractor_adaptive_cell_size.mat