function adj_cells = adjcells(ncell, N)
% Find adjacent cells of a given cell that
% are located in the computational region
%
adj_cells = [];
%
range = {};
z = celltoz(ncell, N);
%
for i = 1:length(N)
    range{i} = [z(i)-1, z(i), z(i)+1];
end
%
comb = allcomb(range{:}); % cell coordinates for all combinations
%
[m, ~] = size(comb);
for i = 1:m
    opt = false;
    z = comb(i,:);
    %
    % elinimate illegal cells out of the computational region
    for j = 1:length(N)
        if isempty( find(z(j) == 1:N(j), 1) )
            opt = true;
            break
        end
    end
    %
    if opt
        continue
    else
        adj_cells = [adj_cells; ztocell(z,N')];
    end
end
%
index = (adj_cells == ncell);
adj_cells(index) = []; % rule out the central cell

