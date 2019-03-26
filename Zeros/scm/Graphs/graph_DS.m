% New directed search algorithm for MOP that involves an utopia vector. A
% cell-space jacobian approximation is made to evaluate a proper descent
% direction with the interpretation of cell space. Evolutionary in
% parameter space can be then determined by taking the linear combination
% of bases acquired from the neighbourhood cells (except the taboo cells).
% The image cell is picked whose base vector has the minimum angle with the
% point-wise evolutionary direction. Like our previous modification, this
% code is also written with the memory saving taken into consideration.
% By: Free Xiong: 2013/06/12
%
% ------------------------------------------------------------------------
% Further modiciation is mading to reduce the size of Fe, Xc and SCM in
% order to save considerable memory during high-dimension calculation. All
% exported arrays are with the same length of S. This will dramatically
% save a lot of memory for high-dimension optimization and largely avoids
% the occurance of "out-of-memory" error.
% By: Free Xiong: 2013/05/08
% -------------------------------------------------------------------------
function [S, SCM, Xc] = graph_DS(dysy, ineqcons, ub, lb, N, S)
%
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
    Fe(i, :) = abs(dysy(xa));
    Xc(i, :) = xa;
end
%
% construct gradient-free ant-colony inspired SCM
for i = 1:length(S)
    cs = S(i);
    neighbours = adjcells(cs, N);
    %
    % judge the validation of cs and its neighbour first
    z = celltoz(cs, N);
    x_cs = ztox(z, h, lb);
    if ineq(ineqcons, x_cs, cs) == 0 % taboo cell
        SCM(i) = 0;
        continue
    else
        % rule out taboo neighbour cells and conduct dominant check to find
        % self-attracted cells
        jRmv = [];
        dominated = false;
        for j = 1:length(neighbours)
            k = find(neighbours(j)==S, 1);
            z = celltoz(neighbours(j), N);
            x = ztox(z, h, lb);
            %
            if isempty(k)
                continue
            else
                if ineq(ineqcons, x, neighbours(j)) == 0 % taboo neighbour
                    jRmv = [jRmv; j];
                    continue
                else
                    if all( Fe(k,:) <= Fe(i,:) )
                        dominated = true; % cs is dominated by its neighbour cell
                    end
                end
            end
        end
        neighbours(jRmv) = [];
        %
        if ~dominated
            SCM(i) = cs; % self-attracted case
        else
            d = -dysy(x_cs); % utopia vector
            [v, B, ~] = Steer_and_Map(x_cs, neighbours, d, dysy, N, lb, ub);
            %
            % image cell will be taken as the one with closest directional
            % angle with the evolutionary direction "v"
            [~, index] = max(cos(v'*B));
            SCM(i) = neighbours(index);
        end
    end
end