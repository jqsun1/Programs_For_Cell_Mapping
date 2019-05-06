% Ant colony inspired mapping construction
% one cell still has several of its adjacent
% cells as potential image cells, but only
% the one with biggest probability will win
%
% ------------------------------------------------------------------------
% Further modiciation is mading to reduce the size of Fe, Xc and SCM in
% order to save considerable memory during high-dimension calculation. All
% exported arrays are with the same length of S. This will dramatically
% save a lot of memory for high-dimension optimization and largely avoids
% the occurance of "out-of-memory" error.
% By: Free Xiong: 2013/05/08
% -------------------------------------------------------------------------
function [S, SCM, Xc, Fe] = graph_free_antcolony(dysy, ineqcons, ub, lb, N, S)
    % dysy here for zero finding problems is the I/O file. In this file, no
    % gradient information which forms the underlying dynamics appears.
    
    h = (ub - lb) ./ N;                 %Computes size of cell in ith
    nz = prod(N);                       %Computes number of cells
    if nargin < 6
        S = (1:nz);
    end
    
    Fe = zeros(length(S),length(dysy(ub)));
    Xc = zeros(length(S),length(N));
    SCM = zeros(length(S),1);
    %
    % evaluate function values of each cell in S
    for i = 1:length(S)
        cs = S(i);
        z = celltoz(cs, N); % Compute coordinates of cs        
        xa = ztox(z,h,lb);
        Fe(i, :) = dysy(xa);
        Xc(i, :) = xa;
    end
    %
    % construct gradient-free ant-colony inspired SCM
    for i = 1:length(S)
        cs = S(i);
        %
        neighbours = adjcells(cs, N);
        %
        dominated = false;
        dFe = []; % function value differences norm
        C = []; % potential image cells
        %
        for j = 1:length(neighbours)
            k = find(neighbours(j)==S, 1);
            %
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
            % different image cell selection schemes are proposed here
            % -------------------------------------------------------
%             [~, index] = min(P);
            [~, index] = max(P);
            dest_cell = C(index); % image cell pick up strategy
            %
%             dest_cell = Roulette_wheel(C,P);
            % -------------------------------------------------------
            %
            z = celltoz(dest_cell,N);
            x = ztox(z, h, lb);
            if ineq(ineqcons, x, dest_cell) == 0 % dest_cell is a taboo cell
                dest_cell = 0;
            end
            SCM(i) = dest_cell;
        else
            z = celltoz(cs, N);
            x = ztox(z, h, lb);
            if ineq(ineqcons, x, cs) == 0 % taboo cell
                SCM(i) = 0;
            else
                SCM(i) = cs;
            end
        end
        %
    end
end