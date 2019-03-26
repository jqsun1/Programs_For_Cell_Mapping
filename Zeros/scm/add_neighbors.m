function S_expand = add_neighbors(S, N, lb, ub)
% -------------------------------------------------------------------------
% Calculate the neigbhbor cells of each cell in set S and expand the set
% with these neighbor cells
% -------------------------------------------------------------------------
h = (ub - lb)./N;
S_expand = S;
for i = 1: length(S)
    cs = S(i);
    neigh = adjcells(cs, N);
    for j = 1: length(neigh)
        z = celltoz(neigh(j), N);
        x = ztox(z, h, lb);
        %
        if ~ismember(neigh(j),S_expand)
            S_expand = [S_expand; neigh(j)];
        end
    end
end
S_expand = unique(S_expand);