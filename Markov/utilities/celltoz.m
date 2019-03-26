function z = celltoz(cell, N)
    cell = cell - 1;
    for i=1:length(N)
        coord(i) = rem(cell, N(i)) + 1;
        cell = fix(cell/N(i));        
    end
    z = coord';
end