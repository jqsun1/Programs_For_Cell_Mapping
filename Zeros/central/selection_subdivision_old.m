function [S_new, N_new, counter] = selection_subdivision_old(f, iter, S_old, N_old, Nm, ...
    div, dsys, g, ub, lb, selection_mode, Q_old)
% -------------------------------------------------------------------------
% Impliment subdivision and selection operations. Also, we have different
% modes on the selection scheme.
%
% Input argument (selected):
%         Q_old:   Block numbers of S_old cell set, Q_old can be empty when
%                  this code is used to conduct the first refinement.
% -------------------------------------------------------------------------
%
terminate = false;
counter = 0;
N1 = N_old;
%
% cell space that has been classified (1st refined set)
S_init = S_old;
N_init = N_old;
%
while true
    %
    % subdivision and filtering
    S_old = reshape(S_old,length(S_old),1);
    [S_new, N_new] = refine(S_old, lb, ub , N1, div); % subdivision from S_old
%     S_new = MOP_examine(f, lb, ub, N_new, S_new); % filtered S_new in terms of MOP
    %
    % selection
    switch selection_mode
        case 1
            % forward shoot Nm times and calculate the intersection of these sets
            C0 = S_new;
            for i = 1:Nm
                [~, C1, ~] = graph_dyn(dsys, g, ub, lb, N_new, C0);
                C0 = intersect(C1, C0);
            end
            % selection, retain only invariant cells after Nm shooting
            S_new(~ismember(S_new, C0)) = [];
        case 2   
            % sequence selection
            [~, C1, ~] = graph_dyn(dsys, g, ub, lb, N_new, S_new);
            S_new = Forward_shoot_eliminate(S_new, C1, Nm);
        case 3
            % refined scm attractor (from MOP)
            [~, C1, ~] = graph_dyn(dsys, g, ub, lb, N_new, S_new); 
            [Gr2, Pe2, St2] = search(C1, N_new, S_new);
            cells = S_new;
            S_new = S_new(Gr2~=1 & St2==0);
%             S_new = add_neighbors(S_new, N_old, lb, ub); % add neighbors
        case 4
            % Directed search via 'min|f|' MOP
            fmop = @(x)abs(f(x));
            [~, SCM, ~, ~] = graph_free_antcolony(fmop, [], ub, lb, N_new, S_new);
            [Gr, ~, St] = search(SCM, N_new, S_new);
            S_new = S_new(Gr~=1 & St==0);
    end
    %
    % check block property of S_new (only for successive iterations)
    if ~isempty(Q_old)
        Q_new = block_check(S_new, N_new, ub, lb, S_init, N_init, Q_old);
        Q_new = unique(Q_new);
        if length(Q_new) ~= length(unique(Q_old))
            terminate = true;
        end
    end
    %
    if terminate
        % if miss block, use the previous results
        S_new = S_old;
        N_new = N1;
        break
    end
    %
    % iteration update
    counter = counter + 1;
    if counter == iter
        break
    else
        S_old = S_new;
        N1 = N_new;
    end
end

% keyboard