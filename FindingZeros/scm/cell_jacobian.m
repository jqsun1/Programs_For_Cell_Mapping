function J = cell_jacobian(f, x, h)
% -------------------------------------------------------------------------
% Jacobian matrix in cell space. The matrix is constructed by doing small
% perturbation with different cell space revolution.
% -------------------------------------------------------------------------
% h = (ub - lb)./N;
a = 0.01; % perturbation parameter
X = repmat(x,1,length(h)) + a*diag(h); % x must be input as a column vector
fx = f(x);
%
J = zeros(length(f(h)),length(h));
for i = 1:size(J,1)
    for j = 1:size(J,2)
        fy = f(X(:,j));
        J(i,j) = (fy(i)-fx(i))/a/h(j);
    end
end