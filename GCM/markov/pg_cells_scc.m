function [pg, gr, type, g] = pg_cells_scc(P,S)
% -------------------------------------------------------------------------
% Using strongly connected components and closed communicating classes to
% identify persistent groups. Tarjan's SCC decomposion and DFS are the main
% ingredients of the searching algorithm. Later modification by integrating
% DFS into Tarjan can hopefully reduce the complexity to O(E+V). Note the
% PG here follows the definition of forward closed only.
%
% Input arguments:
%       P:    Markov chain or graph in sparse matrix
%       S:    Cell set (default is 1:size(P,1))
%
% Output argument:
%      pg:    Invariant set indicator, both stable and unstalbe
%      gr:    Group number of an invariant set, both stable and unstalbe
%    type:    type-1, stable solutions; type-2 unstable solutions
%      g:     Total number of groups, both stable and unstalbe
%
% Note: For the extraction of stable solutions, use the line like this:
%       find(pg==1 & type==1)
%       unique(gr(pg==1 & type==1)) to see the group number of stable sol
%
% Note: This code also illustrates the difference of persistent groups
% acquired from the presented approach and the backward searching.
%
% By: Free Xiong; 2014-08-14, modified on 2014-09-19
% -------------------------------------------------------------------------
if nargin < 2
    S = (1:size(P,1))';
end
[scc_num, scc] = graphconncomp(P);
pg = zeros(length(S),1);
gr = zeros(length(S),1);
type = zeros(length(S),1);
g = 1;
nodes = 1:size(P,1);
for i = 1:scc_num
    scci = nodes(scc==i);
    if(length(scci)==1)
        imgs = find(P(scci,:)~=0);
        if(length(imgs)==1 && imgs==scci)
            % absorbing cell, must be stable PG
            disp(['Absorbing cell found! Assign group number: ',num2str(g)])
            pg(scci) = 1;
            gr(scci) = g;
            type(scci) = 1;
            g = g+1;
        elseif(length(imgs)~=1 && ismember(scci,imgs))
            % single saddle cell, unstalbe
            disp(['Single saddle cell found! Assign group number: ',num2str(g)])
            pg(scci) = 1;
            gr(scci) = g;
            type(scci) = 2;
            g = g+1;
        end
    else
        to = graphtraverse(P,scci(1));
        if length(to)==length(scci)
            % closed communicating class with forward closed only (believed to be PGs)
            disp(['Stable invariant set found! Assign group number: ',num2str(g)])
            pg(scci) = 1;
            gr(scci) = g;
            type(scci) = 1;
            g = g+1;
        else
            % open communicating classes, unstalbe solutions
            disp(['Unstable solution found! Assign group number: ',num2str(g)])
            pg(scci) = 1;
            gr(scci) = g;
            type(scci) = 2;
            g = g+1;
        end
    end
end
g = g-1; % total number of pg (stable and unstable)