clear all
close all
clc
% -------------------------------------------------------------------------
% Test the gcm/scm hybrid approach for covering set finding. The idea is to
% locate several starting cells by using selected scm from gcm. Then,
% traverse the gcm-digraph from these cells as starting cell and calculate
% the sets those cells lead to and are being led. The expected covering set
% is the union of all those sets.
% -------------------------------------------------------------------------
wq = cputime;
%
% up to now the most promising combination is (dopt = 1, sopt = 1 or 7)
sopt = 7; % stepsize control (stepsize.m)
dopt = 1; % local iteration control (dynamic_system.m)
opt = false; % normalization control

mapping_mode = 2; % 1, scm; 2, selected scm; 3, gcm
selection_mode = 3; % 1, set intersection; 2, sequence selection; 3, refined scm

max_iter = 5; % max subdivision times
Nm = 20; % forward shooting times

N_old = [50; 50];
div = 2*ones(length(N_old),1); % (2n+1)

% test on duffing example (3 or 6)
[lb, ub, dsys, g, f] = ProblemDef(3, N_old, dopt, sopt, opt); % problem definition with discrete dyn
GCM = gcmDataBase_aug(N_old, lb, ub, div, dsys);
DG = gcm2graph(GCM); % directed graph
SCM = SelectedSCM(N_old, GCM);
[Gr, Pe, St] = search(SCM, N_old);
cells = find(Gr~=1 & St==0); % these cells must be part of cyclic structures

% traverse the gcm graph and reveal all the cyclic structures
cover = [];
for startcell = cells
    from = graphtraverse(DG, startcell);
    to = graphtraverse(DG', startcell); % inverse mapping
    loop = intersect(from, to);
    loop = union(loop,startcell);
    cover = union(cover,loop); % expand the covering set
end
cover = unique(cover);

time = cputime - wq
PlotCells(cover, lb, ub, N_old, 'red', 'EdgeOn')