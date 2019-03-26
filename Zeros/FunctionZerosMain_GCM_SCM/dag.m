function [G, S, C] = dag(G)
% -------------------------------------------------------------------------
% Compact the nonDAG to a DAG, where the new graph has each node as the
% subgraph No. of the orginal graph.
% -------------------------------------------------------------------------
Nc = size(G,1);
if ~graphisdag(G)
    % nonDAG to DAG
    [S, C] = graphconncomp(G);
    cs = 1:Nc;
    row = [];
    col = [];
    for i = 1:S
        cells = cs(C==i); % cells in the 'ith' subgraph
        for j = 1:length(cells)
            ncell = cells(j);
            for k = 1:Nc
                if G(ncell,k) && ismember(k,cs(C~=i))
                    row = [row; i];
                    col = [col; C(k)];
                end
            end
        end
    end
    %
    dim = unique([row, col],'rows');
    G = sparse(dim(:,1),dim(:,2),true,S,S);
    if ~graphisdag(G)
        error('Compacted graph is still not a DAG!')
    end
end