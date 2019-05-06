function GCM = DS_gcmDatabase(f, ineqcons, ub, lb, N, S)
% -------------------------------------------------------------------------
% Construct GCM based MOP gradient-free mapping. This code has the memory
% save function that is compatiable with all other SCM or GCM storages used
% in previous code sets. This code is actually the successor the SCM
% direct search file 'graph_free_antcolony.m'. Unlike that direct searching
% algorithm, where the maximum function difference norm is taken as the
% unique image cell selection criteria, the GCM based mapping discard that
% criteria and accepts all image cells with all-lower function values.
% Correspondingly, a descent rate, which is a probability like value, will
% be assigned to each image cell. Though this is not used in this code set,
% we set this variable to meet the data structure in consistency with all
% other GCM-associated codes I wrote.
% -------------------------------------------------------------------------
h = (ub - lb) ./ N;                 %Computes size of cell in ith
nz = prod(N);                       %Computes number of cells
if nargin < 6
    S = (1:nz);
end

Fe = zeros(length(S),length(f(ub)));
Xc = zeros(length(S),length(N));
GCM = cell(length(S),3);
%
% evaluate function values of each cell in S
for i = 1:length(S)
    cs = S(i);
    z = celltoz(cs, N); % Compute coordinates of cs
    xa = ztox(z,h,lb);
    Fe(i, :) = f(xa);
    Xc(i, :) = xa;
end
%
% construct gradient-free GCM cell array
for i = 1:length(S)
    cs = S(i);
    neighbours = adjcells(cs, N);
    dominated = false;
    dFe = []; % function value differences norm
    C = []; % potential image cells
    %
    for j = 1:length(neighbours)
        k = find(neighbours(j)==S, 1);
        if isempty(k)
            continue
        else
            if all( Fe(k,:)<=Fe(i,:) )
                dominated = true; % cs is dominated by its neighbour cell
                dFe = [dFe; norm(Fe(k,:)-Fe(i,:))];
                C = [C; neighbours(j)];
            end
        end
    end
    %
    if dominated
        P = dFe/sum(dFe);
        %
        % check inequality satisification for each GCM image cell
        for j = 1:length(C)
            z = celltoz(C(j),N);
            x = ztox(z, h, lb);
            C(j) = ineq(ineqcons, x, C(j));
        end
        % construct GCM in cell array form
        GCM{i,1} = length(C);
        GCM{i,2} = C;
        GCM{i,3} = P;
    else
        z = celltoz(cs, N);
        x = ztox(z, h, lb);
        cs = ineq(ineqcons, x, cs);
        GCM{i,1} = 1;
        GCM{i,2} = cs;
        GCM{i,3} = 1;
    end
end