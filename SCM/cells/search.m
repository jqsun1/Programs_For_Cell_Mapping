function [Gr, Pe, St] = search(C, N, iter)
% New search programm implements unravelling algorithm wirtten by Free
%
% -------------------------------------------------------------------------
% Further modification reducing memory cost is carried out by reducing the
% size of all imported data with the same as "iter". Those cells beyonded
% the set of refined region will be abandoned and are not accounted in the
% memory allocation. Likewise, Gr, Pe and St arrays are reduced by only
% concerning the cells within the refined region. An indexing array is
% introduced to store these feasbile cells.
%
% By: Free Xiong; 2013/05/08
% -------------------------------------------------------------------------
%
nz = prod(N);
if nargin < 3
    ploop = nz;
    iter = 1:nz+1;          % Sink cell is the last entry of cell array, always! JQS
else
    ploop = length(iter);
    iter = [iter,nz+1];     % Add sink cell to the end of a partial array  JQS
end
%
Gr = zeros(ploop,1); % group number
Pe = zeros(ploop,1); % period number
St = zeros(ploop,1); % step number
% Initialize the sink cell
g = 1; % initial group No. for sink cell,
Gr(ploop+1) = g;      % group number
Pe(ploop+1) = 1;      % period number
St(ploop+1) = 0;      % step number
%
% Main frame of unravelling algorithm
for p = 1:ploop
    % loop over cell space "iter"
    cell_start = iter(p);
    if(Gr(cell_start) == 0)
        Gr(cell_start) = -1;
        flag = true; % indicator of each sequence
        SeqIndex = cell_start;
        cell_old = cell_start;
    else
        flag = false;
    end
    %
    while flag
        n = C(cell_old);        % Find image of current cell p
        switch Gr(n) % check image cell status
            case 0 % virgin cell
                Gr(n) = -1; % processing
                cell_old = n; % keep shooting
                SeqIndex = [SeqIndex; n];               
            case -1 %  Find a loop
                % identify j for periodicity
                SeqIndex = [SeqIndex; n];
                j = find(SeqIndex == n);
                j(j==length(SeqIndex)) = [];
                %
                g = g+1; % increase group No.
                Gr(SeqIndex) = g;
                Pe(SeqIndex) = length(SeqIndex)-j;
                St(SeqIndex(j:end)) = 0; % k = j:length(SeqIndex)
                St(SeqIndex(1:1:j-1)) = j-1:-1:1;
                %
                flag = false;
            otherwise % Hit the processed cells
                SeqIndex = [SeqIndex; n];
                Gr(SeqIndex(1:end-1)) = Gr(n);
                Pe(SeqIndex(1:end-1)) = Pe(n);
                St(SeqIndex(1:end-1)) = St(n)+length(SeqIndex)-1: -1: St(n)+1;
                %
                flag = false;
        end
    end
end