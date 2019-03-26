function Q = stack_forw_back(GCM, ie, tgt, S, mode)
% -------------------------------------------------------------------------
% Calculate cell indices of a set that a certain cell leads to or lead to
% that particular cell. All cell-associated information is record in the
% index form of a larger cell set 'S', where cell number rather than their
% indices in that array are assigned.
%
% Input arguments:
%      ie:    particular cell index that starts or ends
%     tgt:    target set, which is supposed to be a subset of 'S'
%      S:     basic space that all cells live in, stored with cell numbers
%     mode:   forward or backward shooting
%
% Output argument:
%      Q:     cell indices that 'ie' leads to or being led, the final
%             results exclude the starting or ending cell index 'ie'
% -------------------------------------------------------------------------
Q_old = []; % indices of cells 'ie' leads to
i_old = ie;
while true
    flag = false;
    for i = i_old
        if strcmp(mode, 'forward')
            img = GCM{i,2};
        elseif strcmp(mode, 'backward')
            img = GCM{i,4};
        else
            error('Either a forward or backward mode should be assigned!')
        end
        %
        if all(ismember(img,tgt)) == 0
            % check if breaking chain exists
            flag = true;
            break
        else
            % store all img/p_img cell indices in 'tgt' set
            i_new = find( ismember(S,intersect(tgt,img)) );
            Q_new = [Q_old; i_new'];
        end
    end
    %
    % emergency break due to chain breaking
    if flag
        break
    end
    %
    % check Q set expansion
    if isempty(setdiff(Q_new,Q_old))
        break
    else
        i_old = unique(i_new);
        Q_old = unique(Q_new);
    end
end
Q = unique(Q_old); % indices 'ie' leads to or being led
Q = setdiff(Q,ie);