function V_new = backward_traj2(V_old, C)
% -------------------------------------------------------------------------
% Backward searching of the set V_old by using just one step forward simple
% cell mapping array C. Note this is the one step backward set expansion.
%
% Input arguements:
%      V_old:     Old set remain to be mapped backwardly
%        C:       Simple cell mapping
%
% Output argument:
%      V_new:    One step backwardly mapped set
% -------------------------------------------------------------------------
%
V_new = [];
for i = 1:length(V_old)
    cell = V_old(i);
    pre_img = find(C==cell);
    V_new = [V_new; pre_img];
end
%
V_new = unique(V_new);