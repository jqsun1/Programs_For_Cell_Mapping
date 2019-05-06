function G = scm2graph(C, iter)
% -------------------------------------------------------------------------
% Transfer simple cell mapping to a graph, pre-image cells of sink cells
% are regarded as attractors in the graph. The graph stores the mapping of
% indices rather than cell number.
% -------------------------------------------------------------------------
row = [];
col = [];
for node = 1:length(iter)
    if C(node) ~= 0
        row = [row; node];
        col = [col; find(C(node)==iter)];
    end
end
G = sparse(row,col,true,length(iter),length(iter));