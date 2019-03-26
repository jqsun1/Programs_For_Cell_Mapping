function [iter, V] = attraction_basin(S, N, lb, ub, dyn_scm)
% -------------------------------------------------------------------------
% Find the attraction basin of a given attractor. The attraction basin is
% the influential region that all cells within this set are absorbed into
% the given attractor. Backward searching algorithm is used here to find
% the reachable set in state space. Note this is the second stage of
% nonlinear system analysis, cell space partition and attractor have all
% being found during the previous stage.
%
% Input arguements:
%       S:      Attractor of a nonlinear system
%   N, lb, ub:  Cell space info of the attractor partition
%    dyn_scm:   Simple Cell Mapping dynamical system
%
% Output arguments:
%     iter:     Total number of backward steps
%      V:       Final reachable set, namely, the attraction basin
% -------------------------------------------------------------------------
%
% Find the vicinity of S as the starting set
V = [];
for i = 1:length(S)
    cell = S(i);
    neig = adjcells(cell, N);
    neig(ismember(neig,S)) = [];
    V = [V; neig];
end
V = unique(V);
%
% Main loop over the backward steps: expand set V
flag = true;
V_old = V;
iter = 0;
[~, C, ~] = graph_dyn(dyn_scm, [], ub, lb, N);
%
while flag
%     V_new = backward_traj(V_old, N, lb, ub, dsys, tf);
    V_new = backward_traj2(V_old, C);
    %
    if all(ismember(V_new,V) == 1)
        break
    else
        V = [V; V_new];
        V = unique(V);
        V_old = V_new;
        iter = iter + 1;
    end
end