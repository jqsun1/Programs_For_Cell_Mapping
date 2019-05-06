function [Q, P] = Adaptive_attractor(B0, N_old, lb, ub, dsys, iter, div)
% -------------------------------------------------------------------------
% Use adaptive scheme to construct the strange attractor. The attractor
% allows the existence of cells with different sizes. Each type of cell is
% stored in a struct variable with the information of cell indices and
% their corresponding cell space partition. The code structure is like a
% binary tree which involves the subdivision and selection operation only
% on a handful of cells with necessity to do so. A classification operation
% at each subdivision iteration is incorporated with the help of
% selection. Hopefully, this new approach could save a large amount of
% computational effort.
%
% Input arguments:
%      B0:    original set covers the attractor (PG set with large cells)
% N_old, lb, ub:   original cell space info
%     dsys:   point-to-point dynamic system
%     iter:   times for subdivision iteration
%      div:   subdivision during each time of a cell (2n+1)
%
% Output argument:
%      Q:   Final attractor with cells of different sizes. Q is stored in an
%           array with each element a struct variable. Cell indices and
%           corresponding cell space partition info is included for further
%           post processing.
%      P:   The complement set of S(i) at each time. P is also stored in
%           the similar way as S and will be used to test the convergence
%           of the algorithm
% -------------------------------------------------------------------------
Q = [];
P = [];
counter = 0;
flag = true;
%
B = B0; % PG set need to be refined all
while flag
    % subdivision on set B, which is classified as discarded category
    [S, N_new] = refine(B, lb, ub, N_old, div);
    h = (ub - lb)./N_new;
    %
    % classification: determine which cell need to be further refined
    % note that S == union(B,R)
    B = []; % cells need to be refined
    R = []; % cells retained as the current size
    %
    if isempty(Q)
        [~, C2, ~] = graph(dsys, ub, lb, N_new, S');
        %
        % classification over the first refinment
        B = S(~ismember(S, C2));
        R = S(ismember(S, C2));
    else
        for i = 1: length(S)
            % return to point-wise judgement
            cell = S(i);
            z = celltoz(cell, N_new);
            x1 = ztox(z, h, lb);
            x2 = dsys(x1);
            %
            % see whether point x2 falls into the evolving invariant set
            %
            % we already have cells with various sizes
            indicator = false;
            %
            for j = 1: length(Q)
                N = Q(j).division;
                members = Q(j).members;
                %
                z = xtoz(x2, h, lb);
                cell = ztocell(z, N);
                %
                if ismember(cell, members)
                    R = [R; S(i)]; % retained at this level of cell
                    indicator = true;
                    break
                else
                    continue
                end
            end
            %
            if ~indicator
                % B set is expanded if the above loop has been completely
                % executed, indicator here is used to detect the execution
                % complete degree.
                B = [B; S(i)];
            end
        end
    end
    %
    % record information about the classification
    Retain.members = R;
    Retain.division = N_new;
    Q = [Q; Retain];
    %
    Refine.members = B;
    Refine.division = N_new;
    P = [P; Refine];
    %
    % update
    counter = counter + 1;
    N_old = N_new;
    %
    % stopping criteria (other criterions can be added here)
    if counter == iter
        flag = false;
    end
end