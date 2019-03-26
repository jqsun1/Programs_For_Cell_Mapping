function [x, lambda] = Stability_check(S_new, N_new, lb, ub, f, dsys)
% -------------------------------------------------------------------------
% At the final stage of computing, check statbility of acquired isolated
% solutions. First, exact solution in point-wise space is solved near each
% isolated block. Then, the Jacobian matrix at these points are evaluated
% to test their stability. Note this post-processing is only suitable for
% scattered solutions. For problems with clustered solutions, it will be
% very costly to run this program
%
% Input arguments:
%     S_new, N_new:  final cell solutions and the partition
%        lb, ub:     boundaries
%      S1, N1, Q1:   all blocks within S1 (original coarse space with N1 partition)
%          f:        nonlinear functions
%        dsys:       dynamic system of iterative scheme
%
% Output arguments:
%         x:         point-wise solutions in each block
%       lambda:      eigenvalue of each jacobian at x
% -------------------------------------------------------------------------
%
Q_new = block_classification(S_new, N_new); % classify new blocks in N_new
Q = unique(Q_new);
h = (ub - lb)./N_new;
%
x = zeros(length(Q),length(h));
lambda = zeros(length(Q),1);
%
for i = 1:length(Q)
    cells = S_new(Q_new==i); % cells in S_new belongs to ith block
    cell = cells(1);
    z = celltoz(cell, N_new);
    x0 = ztox(z, h, lb); % initial guessing
    x(i,:) = fsolve(f,x0); % acquire accurate point-wise solution
    %
    % jacobian matrix in cell space
    J = cell_jacobian(dsys, x(i,:)', h);
%     J = jacobianest(dsys, x(i,:)');
    %
    lambda(i) = max(abs(eig(J)));
end