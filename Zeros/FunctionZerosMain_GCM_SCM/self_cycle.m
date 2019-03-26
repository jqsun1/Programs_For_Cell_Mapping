function G = self_cycle(GCM, N, S)
% -------------------------------------------------------------------------
% Classify a given cell space to different self-cycling groups. The
% algorithm is described in paper [1] on pg. 5.
%
% [1] Improved generalized cell mapping for global analysis of dynamical
% systems. By: Zou HaiLin & Xu JianXue
%
% Input arguments:
%     GCM:   generalized cell mapping database in cell set 'S'. The GCM
%            contains four blocks for each cell, namely, {No. of img cells}
%            {img cells}, {probabilities of img cells}, {pre-image cells}
%      N:    cell space partition
%      S:    initial target set (default as the whole cell space)
%
% Output arguments:
%      G:    self_cycle groups in cell space, discard cells are assigned
%            with g=-1, which implies they are transient cells
% -------------------------------------------------------------------------
if nargin < 3
    Nc = prod(N);
    S = 1:Nc; % initial target set
    tgt = S;
else
    tgt = S;
end
%
G = zeros(length(S),1);
g = 0;
%
iter = 0;
while ~isempty(tgt)
    iter = iter + 1
    % 1. Conduct cell set delete
    iRmv = [];
    for i = 1:length(tgt)
        img = GCM{S==tgt(i),2}; % note the gcm indices is in accordance with S
        p_img = GCM{S==tgt(i),4};
        if all(ismember(union(img,p_img),tgt)) == 0
            % no pre-image and image cells are in the target set
            iRmv = [iRmv; i];
        end
    end
    d_tgt = tgt;
    d_tgt(iRmv) = [];
    %
    % 2. Apply theorm 2 for individual self-cycle finding
    if isempty(d_tgt)
        break
    else
        % 2.1 choose a proper starting cell and record its cell path
        for i = 1:length(d_tgt)
            flag = false;
            p_img = GCM{S==d_tgt(i),4};
            img = GCM{S==d_tgt(i),2};
            %
            if isempty(intersect(p_img,d_tgt)) || isempty(intersect(img,d_tgt))
                continue
            else
                % 2.2 generate a cell path to find a self-cycling strcture
                path = [];
                is = find(S==d_tgt(i)); % starting cell index
                while true
                    path = [path; is]; % indices of S
                    indices = find(ismember(S,GCM{is,2})); % img cells' indices
                    %
                    % in case there is no image cells in d_tgt set
                    if isempty(intersect(S(indices),d_tgt))
                        break
                    end
                    %
                    if all(ismember(indices,path)) == 1
                        ie = path(end); % 'ie' belongs to a self-cycling structure
                        flag = true;
                        break
                    else
                        for j = indices
                            if ~ismember(j,path)
                                in = j; % choose a img cell index never shown in the path
                                break
                            end
                        end
                        is = in; % keep shooting forward
                    end
                end
            end
            %
            if flag
                break
            end
        end
        %
        % 2.3 use the theorm (b) to extract one complete self-cycling from 'ie'
        Q = stack_forw_back(GCM, ie, d_tgt, S, 'forward');
        P = stack_forw_back(GCM, ie, d_tgt, S, 'backward');
        %
        % no more self-cycling structure can be found
        if isempty(Q) || isempty(P)
            break
        end
        i_cycle = union(union(P,Q),ie); % stored in indices of S
        %
        % 2.4 locate one complete self-cycling and delete related cells
        %     from d_tgt set, meanwhile, update tgt set
        g = g + 1;
        G(i_cycle) = g;
        tgt = setdiff(tgt,S(i_cycle)); % delete procssed cells
    end
end