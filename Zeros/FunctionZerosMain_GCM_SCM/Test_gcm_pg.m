clear all
close all
clc
% test on gcm_pg.m for covering set finding
%
% up to now the most promising combination is (dopt = 1, sopt = 1 or 3)
sopt = 3; % stepsize control (stepsize.m)
dopt = 1; % local iteration control (dynamic_system.m)
opt = false; % normalization control
%
wq = cputime;
N = [16; 16];
[lb, ub, dsys, g] = ProblemDef(7, N, dopt, sopt, opt); % problem definition with discrete dyn
%
div = ones(length(N),1); % 2n+1
GCM = gcmDataBase(N, lb, ub, div, dsys);

%% First stage coarse finding of P-K cells
% [~, C, ~] = graph_dyn(dsys, g, ub, lb, N); % scm with underlying dyn sys
% [~, C, ~] = graph_DS(dsys, g, ub, lb, N_old); % cell space steering
% [~, C, ~] = graph_free_antcolony(dsys, g, ub, lb, N_old); % cell space steering
C = SelectedSCM(N, GCM);

[Gr, Pe, St] = search(C, N);
scm_att = find(Gr~=1 & St==0);

%% Expand the P-K cell set into a complete PG-like set
% [S_pg, ~] = ExtraCell(dsys, div, scm_att, N, lb, ub); % include extra cells to reveal complete info
% time1 = cputime - wq

% [S_pg, S_new, G] = gcm_pg(GCM, scm_att, N);
% time2 = cputime - wq

%%
figure(1)
PlotCells(scm_att, lb, ub, N, 'red', 'EdgeOn')
% hold on
% PlotCells(S_new, lb, ub, N, 'green', 'EdgeOn')
% whos S_new