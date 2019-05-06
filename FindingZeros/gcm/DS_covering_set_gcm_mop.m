function [S, GCM] = DS_covering_set_gcm_mop(f, N, lb, ub, ineqcons, iter, stb)
% -------------------------------------------------------------------------
% Find the covering set using scm/gcm hybrid algorithm. This approach
% starts by picking those P-K cells under the selected scm from gcm and
% find all cyclic structures with those P-K cells involved in. Note for a
% given incomplete sub cell set, the nodes is represented by corresponding
% indices in the cell array.
% -------------------------------------------------------------------------
if nargin < 6
    nc = prod(N);
    iter = 1:nc;
    stb = 'stable';
end
%
% ---------------------------------------------------
% This is the only difference with other gcm/scm
% covering cell finding programs.
GCM = DS_gcmDatabase(f, ineqcons, ub, lb, N, iter);
% ---------------------------------------------------
%
[DG, P] = gcm2graph(GCM, iter); % directed graph
if strcmp(stb, 'unstable')
    DG = DG';
    P = P';
end
%
if graphisdag(DG)
    error('Current set is acyclic!');
else
    SCM = SelectedSCM_graph(N, DG, P, iter, 'max');
    [Gr, Pe, St] = search(SCM, N, iter);
    cell_inds = find(Gr~=1 & St==0); % these cells must be part of cyclic structures
    %
    if isempty(cell_inds)
        % in case there is no P-K cells being captured, traverse whole set
        warning('No P-K cells are found under SCM, traverse the whole set')
        cell_inds = 1:length(iter);
    end
    %
    % traverse the gcm graph and reveal all the cyclic structures
    cover_inds = [];
    for startnode = cell_inds
        from = graphtraverse(DG, startnode);
        to = graphtraverse(DG', startnode); % inverse mapping
        loop = intersect(from, to);
        loop = union(loop,startnode);
        cover_inds = union(cover_inds,loop); % expand the covering set
    end
    cover_inds = unique(cover_inds);
    S = iter(cover_inds); % output cell numbers
end