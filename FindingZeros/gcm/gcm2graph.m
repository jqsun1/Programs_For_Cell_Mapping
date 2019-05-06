function [G, P] = gcm2graph(GCM, iter)
% -------------------------------------------------------------------------
% Transfer GCM database to a matlab recognized graph representation. Sparse
% matrix over the cell space is used to generate the graph with a logical
% bit representing the location of image cell. Corresponding transitional
% probability from a domain cell to its image cells are also stored in
% another sparse matrix.
%
% Output arguments:
%       G:    sparse matrix of image cells' location (logical)
%       P:    sparse matrix of transitional probabilities (real)
% -------------------------------------------------------------------------
Nc = size(GCM,1);
if nargin < 2
    iter = 1:Nc;
end
%
rows = [];
cols = [];
Pro = [];
for i = 1:length(iter) % iter stores all cell numbers
    I = GCM{i,1};
    for j = 1:I
        % we don't take any edges that lead to sink cells
        if GCM{i,2}(j) ~= 0
            rows = [rows; i];
            cols = [cols; find(GCM{i,2}(j)==iter)];
            Pro = [Pro; GCM{i,3}(j)];
        end
    end
end
G = sparse(rows,cols,true,length(iter),length(iter));
P = sparse(rows,cols,Pro,length(iter),length(iter));