function Q_new = block_check(S_new, N_new, ub, lb, S_old, N_old, Q_old)
% -------------------------------------------------------------------------
% Designate block number of a newly refined set which is based on the
% coarse cell set. The code is used to check whether an iteration has
% killed potential solutions in a too refined cell set.
%
% Input arguments:
%     S_new, N_new:   Newly refined cell set and its cell space partition
%        ub, lb:      Upper and lower bound of computational region
%     S_old, N_old:   Old coarse cell set that has different blocks
%         Q_old:      Block numbers of each cell in S_old
%
% Output argument:
%         Q_new:      Block number of each cell in S_new
% -------------------------------------------------------------------------
%
% Q_new = zeros(length(S_new),1);
cs = zeros(length(S_new),1);
h_new = (ub - lb)./N_new;
h_old = (ub - lb)./N_old;
%
for i = 1:length(S_new)
    % translate to coarse grid
    cell = S_new(i);
    z = celltoz(cell, N_new);
    x = ztox(z, h_new, lb);
    z = xtoz(x, h_old, lb);
    cs(i) = ztocell(z, N_old);
end
%
cs = unique(cs);
Q_new = Q_old(ismember(S_old,cs));