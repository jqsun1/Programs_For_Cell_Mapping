function [S_new, N_new, counter, mapping] = selection_subdivision(iter, S_old, N_old, Nm, ...
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
            % forward shoot Nm times and calculate the intersection of these sets
            C0 = S_new;
            for i = 1:Nm
                [~, SCM, ~] = graph(dsys, ub, lb, N_new, C0);
                C0 = intersect(SCM, C0);
            end
            % selection, retain only invariant cells after Nm shooting
            S_new(~ismember(S_new, C0)) = [];
        case 2   
            % repeat the process of GCM-SCM hybrid to find smaller cycles
%             S_new = add_neighbors(S_new, N_new, lb, ub); % add neighbors
            [S_new, GCM] = covering_set_gcm(N_new, lb, ub, div, dsys, S_new, stb);
        case 3
            % discard those with image cells all out
            GCM = gcmDataBase_aug(N_new, lb, ub, div, dsys, S_new);
            iRmv = [];
            for i = 1:size(GCM,1)
                img = GCM{i,2};
                if isempty(intersect(img,S_new))
                    iRmv = union(iRmv,i);
                end
            end
            S_new(iRmv) = [];
        case 4 % This is the most reliable selection algorithm
            % discard those with pre-image cells all out (gcm)
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
        case 5
            % discard those with pre-image cells all out (scm)
            [~, SCM, ~] = graph(dsys, ub, lb, N_new, S_new);
            DG = scm2graph(SCM, S_new);
            iDG = DG';
            iRmv = [];
            for i = 1:length(SCM)
                if isempty(S_new(iDG(i,:)))
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