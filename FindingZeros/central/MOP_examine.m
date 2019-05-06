function S_filtered = MOP_examine(f, lb, ub, N, S)
% -------------------------------------------------------------------------
% Examine the neighborhood cells of each cell in a set by checking the
% absoluate values of the function vector. The code aims at reducing the
% number of cells in a set that remain to be processed by SCM unravelling
% algorithm for zero finding, especially in high dimensional space. The
% principle of MOP here is to compare the adjacent cells' function values
% with the central cell. If all the function values are further away from
% zeros, the neighborhood cell might be a transient cell under SCM, and
% there is normally no need to calculate that cell in SCM programs.
% However, the "worse" cell equals transient asertion has no strict proof
% yet. Be careful with the filtered cell space, it might miss potential
% solutions!
% -------------------------------------------------------------------------
h = (ub - lb)./N;
iRmv = [];
for i = 1:length(S)
    cell = S(i);
    z = celltoz(cell, N);
    xc = ztox(z, h, lb);
    %
    neig = adjcells(cell, N);
    neig = neig(ismember(neig,S)); % all neighbor cells examines are in S
    %
    if ~isempty(neig)
        for j = 1:length(neig)
            cs = neig(j);
            z = celltoz(cs, N);
            xj = ztox(z, h, lb);
            % xj is "worse" than xc in terms of zero finding
            if all(abs(f(xj)) > abs(f(xc)))
                iRmv = [iRmv; find(cs==S)];
            end
        end
    end
end
S(iRmv) = [];
S_filtered = S;