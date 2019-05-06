function Q = block_classification(S, N)
% -------------------------------------------------------------------------
% Classify clustered cells into one compact block. This will be used later
% for cell space trust region identification. The classicfication is done
% by a continuation like method, which starts with a unassigned cell and
% shoot forward by finding its neighbourhood cells and repeat the searching
% until a block set is full.
%
% Input arguments:
%         S:    cell set remain to be classified
%         N:    cell space partition
%
% Output argument:
%         Q:    block number array, same length as S
% -------------------------------------------------------------------------
Q = zeros(length(S),1);
blk = 0;
%
while ismember(0, Q)
    q = S(Q==0);
    q_old = q(1); % starting cell of a new block
    blk = blk + 1;
    %
    while true
        % use the continuation like idea to shoot forward
        temp = [];
        for i = 1:length(q_old)
            cell = q_old(i);
            neigh = adjcells(cell, N);
            temp = [temp; neigh];  % all neighbours
        end
        temp = temp(ismember(temp, S));
        q_old_neigh = setdiff(temp, q_old);
        %
        % termination condition
        q_new = [q_old; q_old_neigh]; % expand the current block
        Q(ismember(S, q_new)) = blk;
        if isempty(setdiff(q_new,q_old))
            break
        else
            % data update
            q_old = q_old_neigh;
            q_old = q_new;
        end
    end
end