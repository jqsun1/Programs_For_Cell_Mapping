function [pg, gr, g] = pg_cells_scc2(P,S)
% -------------------------------------------------------------------------
% Using strongly connected components and closed communicating classes to
% identify persistent groups. Note this version only uses one Tarjan's
% algorithm. Persistent groups are those strongly connected components with
% zero out degree in terms of that group. Here, the PGs are treated as the
% SCCs with forward and backward closed.
%
% Input arguments:
%       P:    Markov chain or graph in sparse matrix
%       S:    Cell set (default is 1:size(P,1))
%
% Output argument:
%      pg:    Closed communicating classes, namely, persistent groups
%
% By: Free Xiong; 2014-08-14
% -------------------------------------------------------------------------
if nargin < 2
    S = (1:size(P,1))';
end
[scc_num, scc] = graphconncomp(P);
pg = zeros(length(S),1);
gr = zeros(length(S),1);
g = 1;
nodes = 1:size(P,1);
for i = 1:scc_num
    scci = nodes(scc==i);
    if(length(scci)==1)
        if(P(scci,scci)==1)
            % absorbing cell, must be PG
            pg(scci) = 1;
            gr(scci) = g;
            g = g+1;
        end
    else
        to = graphtraverse(P,scci(1));
        from = graphtraverse(P',scci(1));
        c = intersect(from,to);
        if length(c)==length(scci)
            % closed communicating class (believed to be PGs)
            pg(scci) = 1;
            gr(scci) = g;
            g = g+1;
        end
    end
end
g = g-1; % total number of pg in the gcm