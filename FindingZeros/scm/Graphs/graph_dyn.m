% Output:
% S: cell index array (1-Total No. Cell)
% C: cell-to-cell mapping
% Xc: central point of each in the space
% Fe: f(Xc) of each cell in the space
%
% ------------------------------------------------------------------------
% Further modiciation is mading to reduce the size of Fe, Xc and C in
% order to save considerable memory during high-dimension calculation. All
% exported arrays are with the same length of S. This will dramatically
% save a lot of memory for high-dimension optimization and largely avoids
% the occurance of "out-of-memory" error.
% By: Free Xiong: 2013/05/08
% -------------------------------------------------------------------------
function [S, C, Xc] = graph_dyn(dysy, inequcons, ub, lb, N, S)
    h = (ub - lb) ./ N;                 %Computes size of cell in ith
    nz = prod(N);                       %Computes number of cells
    if nargin < 6
        S = (1:nz);
    end
    C = zeros(length(S), 1);
%     Fe = zeros(length(S),length(f(ub)));
    Xc = zeros(length(S),length(N));
    for i = 1:length(S)
        cs = S(i);
        
        z = celltoz(cs, N); % Compute coordinates of cs        
        xa = lb + h .* z - 1/2 * h;           %Compute center point
        xd = dysy(xa);                        % Point evolution
        z = round((xd - lb) ./ h + 1/2);      %Compute cell coordinates
        ncell = ztocell(z, N);

%         Fe(i, :) = fa;
        Xc(i, :) = xa;
        if ~isempty(find(S == ncell, 1)) %&& (ineq(ineqcons, xd', cs) ~= 0)
            C(i) = ncell; % cs is not a sink cell
        else
            C(i) = 0;
        end                
    end
end