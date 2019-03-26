function P = Markov(dsys, N, lb, ub, S, div)
% -------------------------------------------------------------------------
% Construct the Markov transitional matrix by sampling appraoch. The Markov
% matrix is stored in a spare matrix.
%
% Input Arguments:
%     dsys:    dynamic system for point to point mapping
%  N, lb, ub:  cell space info
%      S:      a given PG-like set for probability estimation
%     div:     sampling subdivision
%
% Output Argument:
%      P:      spare markov transitional matrix
% -------------------------------------------------------------------------
%
% proceed column wisely
rows = [];
cols = [];
s = [];
%
for j = 1:length(S)
    ncell = S(j);
    [I, C, P] = OneCellGCM(ncell, N, lb, ub, div, dsys, S);
    %
    for m = 1:I
        if ismember(C(m),S)
            i = find(C(m)==S);
        else
%             i = length(S) + 1; % for sink cell
            continue
        end
        %
        rows = [rows; i];
        cols = [cols; j];
        s = [s; P(m)];
    end
end
%
% % construct the right lower corner sink cell info
% rows = [rows; length(S)+1];
% cols = [cols; length(S)+1];
% s = [s; 1];
%
% spare matrix storage
% P = sparse(rows, cols, s, length(S)+1, length(S)+1);
P = sparse(rows, cols, s, length(S), length(S));