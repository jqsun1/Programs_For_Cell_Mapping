function unstable_mani_id = unstable_manifold(gr,type,P)
% -------------------------------------------------------------------------
% Locate the unstable manifold from unstable solutions, the idea is to
% traverse the GCM start from the unstable invariant set.
%
% Input arguments:
%       gr:    group number of stable and unstable invariant set
%     type:    stability type, type-1 is stable, type-2 is unstable
%       P:     GCM graph
%
% Output argument:
%  unstable_mani_id: logic array with unstable manifold indicator
%
% By: Free Xiong; 2014-09-24
% -------------------------------------------------------------------------
gr_unstable = unique(gr(gr~=0 & type==2));
unstable_mani_id = zeros(length(P),1);
nodes = [];
for i = 1:length(gr_unstable)
    scci = find(gr==gr_unstable(i));
    disc = graphtraverse(P,scci(1));
    nodes = union(nodes,disc);
end
unstable_mani_id(nodes) = 1;