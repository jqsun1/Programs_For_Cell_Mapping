function S = covering_set(mapping_mode, dsys, g, ub, lb, N_old, div)
% -------------------------------------------------------------------------
% Build mapping according to different requests
% 1. pure scm by recognizing central point
% 2, selected scm from gcm with maximum transition possibility
% 3, gcm persistent group finding
% -------------------------------------------------------------------------
%
switch mapping_mode
    case 1
        % build pure SCM using central points
        [~, C, ~] = graph_dyn(dsys, g, ub, lb, N_old); % scm with underlying dyn sys
        % [~, C, ~] = graph_DS(dsys, g, ub, lb, N_old); % cell space steering
        % [~, C, ~] = graph_free_antcolony(dsys, g, ub, lb, N_old); % cell space steering
    case 2
        % build the selected SCM from GCM
        GCM = gcmDataBase(N_old, lb, ub, div, dsys);
        C = SelectedSCM(N_old, GCM);
    case 3
        % build the GCM directly
        GCM = gcmDataBase(N_old, lb, ub, div, dsys);
        C = SelectedSCM(N_old, GCM);
    case 4
        % max and min probability of selected SCM
        GCM = gcmDataBase(N_old, lb, ub, div, dsys);
        Cmax = SelectedSCM(N_old, GCM, 1:prod(N_old), 'max');
        Cmin = SelectedSCM(N_old, GCM, 1:prod(N_old), 'min');
end
%
% find general covering set based on different mappings
if mapping_mode == 1 || mapping_mode == 2
    % unravelling from SCM, extract invariant set
    [Gr, ~, St] = search(C, N_old);
    scm_att = find(Gr~=1 & St==0);
    S = scm_att;
elseif mapping_mode == 3
    % unravelling from SCM, extract invariant set
    [Gr, ~, St] = search(C, N_old);
    scm_att = find(Gr~=1 & St==0);
    gcm_att = gcm_pg(GCM, scm_att, N_old);
    S = gcm_att;
elseif mapping_mode == 4
    % unravelling from SCM, extract invariant set
    [Gr1, ~, St1] = search(Cmax, N_old);
    scm_att1 = find(Gr1~=1 & St1==0);
    [Gr2, ~, St2] = search(Cmin, N_old);
    scm_att2 = find(Gr2~=1 & St2==0);
    S = [scm_att1; scm_att2];
    S = unique(S);
end