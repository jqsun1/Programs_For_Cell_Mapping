function [S_new, N_new, counter, mapping] = selection_subdivision(iter, S_old, N_old, ...
    div, dsys, ub, lb, selection_mode, stb)
% -------------------------------------------------------------------------
% Impliment subdivision and selection operations. Also, we have different
% modes on the selection scheme.
%
% Input argument (selected):
%         stb:  Find stable (forward) or unstable (backward) attractor
% -------------------------------------------------------------------------
if nargin < 10
    stb = 'stable';
end
counter = 0;
%
while true
    %
    % subdivision
    S_old = reshape(S_old,length(S_old),1); % accept only column vector for refine
    [S_new, N_new] = refine(S_old, lb, ub , N_old, div); % subdivision from S_old
    %
    % selection
    switch selection_mode
        case 1            
            % refined scm attractor (from MOP)
            [~, C1, ~] = graph_dyn(dsys, g, ub, lb, N_new, S_new); 
            [Gr2, ~, St2] = search(C1, N_new, S_new);
            S_new = S_new(Gr2~=1 & St2==0);
        case 2
            % gcm backward searching
            GCM = gcmDataBase_aug(N_new, lb, ub, div, dsys, S_new);
            DG = gcm2graph(GCM, S_new);
            if strcmp(stb, 'unstable')
                DG = DG';
            end
            %
            iRmv = [];
            iDG = DG'; % inversed GCM (either for forward and backward)
            for i = 1:size(GCM,1)
%                 pimg = GCM{i,4};
                pimg = S_new(iDG(i,:)); % use graph to find pre-image
                if isempty(intersect(pimg,S_new))
                    iRmv = union(iRmv,i);
                end
            end
            S_new(iRmv) = [];
    end
    %
    % iteration update
    counter = counter + 1;
    if counter == iter
        break
    else
        S_old = S_new;
        N_old = N_new;
    end
end
%
if ismember(selection_mode,[2,3,4]);
    mapping = GCM;
else
    mapping = SCM;
end