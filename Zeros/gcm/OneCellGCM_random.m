function [I, C, P] = OneCellGCM_random(ncell, N, lb, ub, k, ds, S)
% Construct GCM for a given cell by k-FF sampling method. The random walk
% sampling technique is used to break the curse of dimension somewhat.
%
Nc = prod(N);
if nargin < 7
    S = 1:Nc;
end
%
% x_samples = kFF_sampling(ncell, lb, ub, N, k);
h = (ub - lb)./N;
z = celltoz(ncell,N);
xc = ztox(z, h, lb);
cell_lb = xc - h/2;
n = length(lb);
x_unif = repmat(cell_lb',k,1) + repmat(h',k,1).*rand(k,n); % candidate point set
x_lhs = repmat(cell_lb',k,1) + repmat(h',k,1).*lhsdesign(k,n,'smooth','off'); % latin hypercube sampling
x_samples = [x_unif; x_lhs; xc'];
%
C_original = zeros(size(x_samples,1),1);
for i = 1:size(x_samples,1)
    x_old = x_samples(i,:)';
    x_new = ds(x_old);
    z = xtoz(x_new, h, lb);
    img_cell = ztocell(z, N);
    if ismember(img_cell,S)
        C_original(i) = img_cell;
    end
end
%
C = unique(C_original);
I = length(C);
P = zeros(I,1);
for i = 1:I
    indices = find(C(i)==C_original);
    Mi = length(indices);
    P(i) = Mi/I;
end