function [ipg, igr] = backward_scm(S,C)
% -------------------------------------------------------------------------
% Construct backward scm and sort the closed communicating group from the
% inverse scm. Note generally the inverse scm is a gcm and thus graph
% algorithms are appiled here for analysis.
%
% Input arguments:
%       S:        cell set from 1 to prod(N)
%       C:        simple cell mapping from of S
%
% Output arguments:
%       ipg:      persistent group indicator (0 or 1)
%       igr:      persistent group number of inverse scm, those cells with
%                 ipg == 1 are the PG cells in terms of backward scm.
%
% By: Furui Xiong; 10/28/2015
%--------------------------------------------------------------------------
ipg = zeros(size(S));
igr = zeros(size(S));

% construct the inverse scm and sort out the persistent group, note the
% backward scm is actually a gcm
rows = S;
cols = C;
nc = length(S);
rows(~ismember(cols,1:nc)) = []; % rule out edges to sink cells
cols(~ismember(cols,1:nc)) = [];
G = sparse(rows,cols,true,nc,nc);
P = G'; % inverse scm in graph representation
[scc_num, scc] = graphconncomp(P);

% check the topological property of each strongly connect component
g = 1;
for i = 1:scc_num
    scci = S(scc==i);
    if(length(scci)==1)
        imgs = find(P(scci,:)~=0);
        if(length(imgs)==1 && imgs==scci)
            % absorbing cell, must be stable PG
            ipg(scci) = 1;
            igr(scci) = g;
            g = g+1;
        elseif(length(imgs)~=1 && ismember(scci,imgs))
            % single saddle cell, unstalbe
            ipg(scci) = 1;
            igr(scci) = g;
            g = g+1;
        end
    else
        to = graphtraverse(P,scci(1));
        if length(to)==length(scci)
            % closed communicating class with forward closed only (believed to be PGs)
            ipg(scci) = 1;
            igr(scci) = g;
            g = g+1;
        else
            % open communicating classes, unstalbe solutions
            ipg(scci) = 1;
            igr(scci) = g;
            g = g+1;
        end
    end
end