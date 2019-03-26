function [S, N] = covering_set_gcm_gradual(N, lb, ub, div, dsys, k, maxiter, iter)
% -------------------------------------------------------------------------
% This code the update of its predesour of covering_set_gcm.m, which
% searchs the persitent group under GCM using DFS in graph theory. Unlike
% the previous approach, the GCM covering search is conducted in an
% iterative manner with gradually enhanced resolution. Random sampling kFF
% technique is used to suit higher dimension calculation.
% -------------------------------------------------------------------------
for i = 1:maxiter
    % Refine the incoming set 'iter'
    iter = reshape(iter,length(iter),1);
    if i ~= 1
        [iter, N] = refine(iter, lb, ub, N, div);
    end
    
    % random sampling for GCM
    GCM = gcmDataBase_random(N, lb, ub, k, dsys, iter);
    [DG, P] = gcm2graph(GCM, iter); % directed graph
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
    %
    iter = S;
end