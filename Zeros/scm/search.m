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
    iter = 1:nz;
end
g = 1; % initial group No. for sink cells
%
Gr = zeros(length(iter),1); % group number
Pe = zeros(length(iter),1); % period number
St = zeros(length(iter),1); % step number
%
% Main frame of unravelling algorithm
for p = 1:length(iter)
    % no necessary to extract paths from processed cells !!!
    if Gr(p) ~= 0 && Gr(p) ~= -1
        continue
    end
    % loop over cell space "iter"
    cell_start = iter(p);
    flag = true; % indicator of each sequence
    cell_old = cell_start;
    %
    % loop over a particular sequence with unknown cells
%     SeqIndex = find(cell_start==iter, 1);
    SeqIndex = [];
    while flag
        n = find(cell_old==iter, 1);
        %
        if ~isempty(n)
            switch Gr(n) % check cell status
                case 0 % virgin cell
                    Gr(n) = -1; % processing
                    cell_new = C(n);
                    cell_old = cell_new; % keep shooting
                    SeqIndex = [SeqIndex; n];
                    continue
                case -1 % limit cycle
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
                otherwise % two sequences mingle
                    SeqIndex = [SeqIndex; n];
                    Gr(SeqIndex(1:end-1)) = Gr(n);
                    Pe(SeqIndex(1:end-1)) = Pe(n);
                    St(SeqIndex(1:end-1)) = St(n)+length(SeqIndex)-1: -1: St(n)+1;                    
                    %
                    flag = false;
            end
        else % encounter sink cell
            Gr(SeqIndex) = 1;
            Pe(SeqIndex) = 1;
            St(SeqIndex) = 0;
            flag = false;
        end 
    end
end