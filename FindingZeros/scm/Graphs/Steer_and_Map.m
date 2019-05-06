function [v, V, F] = Steer_and_Map(x0, cell_neig, d, f, N, lb, ub)
% -------------------------------------------------------------------------
% Calculate best evolutionary vector "df" in domian space by letting it
% stay as close as the utopia vector "d". Use the linear combination
% pattern to construct evolutionary vector "dv" in parameter space. Note
% that "div" is acquired by combining all unit vectors point to the central
% point of each cell stored in "cell_neig"
%
% Input arguments:
%       x0:    A given point in parameter space
%   cell_neig: Selected neighbourhood cells around x0
%       d:     A given utpoia vector in domain space (-f(x0))
%       f:     MOP vector field, namely the input-output function
%  N, lb, ub:  SCM partition parameters
%
% Output arguments:
%       V:   Base vectors of neigbourhood cells (each column)
%       F:   Jacobian matrix in cell space (dominant check)
%       v:   Direction in parameter space
% -------------------------------------------------------------------------
%
% store info among those selected neighbourhood cells
xn = zeros(length(cell_neig),length(N));
fn = zeros( length(cell_neig),length(f(x0)) );
h = (ub - lb)./N;

for i = 1:length(cell_neig)
    cs = cell_neig(i);
    z = celltoz(cs, N);
    xn(i,:) = ztox(z, h, lb);
    fn(i,:) = f(xn(i,:)');
end
%
% construct minimal norm matrix
F = zeros(length(f(x0)),length(cell_neig));
V = zeros(length(x0),length(cell_neig)); % used to construct v
for j = 1:length(cell_neig)
    F(:,j) = (f(xn(j,:)') - f(x0))/norm(xn(j,:)'-x0);
    V(:,j) = (xn(j,:)'-x0)/norm(xn(j,:)'-x0);
end
%
% use lambda to construct evolutionary direction in parameter space
lambda = pinv(F)*d; % minimal norm appraoch
v = V*lambda;