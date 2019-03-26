function k_FF_set = kFF_sampling(ncell, lb, ub, N, k)
% -------------------------------------------------------------------------
% k-Furthest First sampling technique in a given cell. The sampling utlizes
% a random walk to explore a confine region with maximum test points
% distribution. The k-FF will be used to sample high-dimension GCM with
% fixed number of testing points.
%
% Input arguments:
%     ncell:   cell number that defines a bounded region
%    lb,ub,N:  cell space partition
%       k:     k points locate in ncell will be sampled
%
% Output argument:
%    k_FF_set: point set represent a cell for GCM construction, k_FF_set is
%              a (k+1)-by-n matrix with cell center included and each
%              sample point list as a row vector
% -------------------------------------------------------------------------
h = (ub - lb)./N;
z = celltoz(ncell, N);
xc = ztox(z, h, lb);
lb = xc - h/2; % lower bound of a cell
ub = xc + h/2; % upper bound of a cell
n = length(lb);
%
% generate (beta*k) random points in ncell
beta = 20;
k1 = beta*k; % selection of k-FF set is among k1 uniformly distributed points

S_old = repmat(lb',k1,1) + repmat(ub'-lb',k1,1).*rand(k1,n); % represent ncell randomly
C = repmat(lb',k1,1) + repmat(ub'-lb',k1,1).*rand(k1,n); % candidate point set
%
% k-FF (k-Furtherest First) procedure
k_FF_set = []; % which is denoted as N in ref [1]
for j = 1:k
    S_expan = [S_old;k_FF_set];
    %
    % Calculate Eqn (1) in ref [1]
    [rows, ~] = size(S_expan);
    MinDist = [];
    for p = 1:length(C)
        % distances between one randomly sampled point and all
        % other points in S_expan
        dists = sqrt( sum((repmat(C(p,:),rows,1) - S_expan).^2,2) );
        MinDist = [MinDist;min(dists)];
    end
    %
    % Find the canditate for k_FF_set from C
    [~, index] = max(MinDist);
    k_FF_set = [k_FF_set; C(index,:)];
    C(index,:) = [];
end
%
% cell center is added for inclusion
k_FF_set = [k_FF_set; xc'];