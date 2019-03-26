function [py, status] = PoincareMap(y, index, c, f, tf, lb, ub)
% -------------------------------------------------------------------------
% Poincare iteration construction as point-to-point mapping. The iteration
% starts from an initial point locates on the Poincare intersection defined
% as {x_index=c}, which means the x_index of the original dynamic system is
% freezed. After a relatively long time integration, the results could
% either relocated on the intersection or not. For either cases, record
% either the re-intersection point or the final point. This point-to-point
% mapping would be furtherly used to construct the cell-to-cell mapping in
% this intersection plane to look at the global attractor for this profile.
%
% Input arguements:
%       y:    Starting point on the {x_index=c} intersection
%     index:  Index No. of the solution vector x_index
%       c:    Intersection position along x_index direction
%       f:    Dynamic system dx/dt = f(x) with function handle
%     lb, ub: Lower and upper bound of Poincare intersection
%      tf:    Integration time of the dynamical system (full dimension),
%             this integration time is ususally relatively long
%
% Output arguement:
%      py:    Poincare iteration point or the final integration point
%    status:  Wether the other Poincare intersection point has be found
% -------------------------------------------------------------------------
N = length(y) + 1;
x0 = zeros(N,1); % N dim state space
%
% extend the initial condition to suit the dynamic system's integration
x0(index) = c;
if index == 1
    x0(index+1:end) = y;
elseif index == N
    x0(1:index-1) = y;
else
    x0(1:index-1) = y(1:index-1);
    x0(index+1:end) = y(index:end);
end
%
% integrate the dynamic system with relatively long period of time
[t, x] = ode45(f,[0 tf],x0); % record all info
% [t, x] = ode45solve(f, 0:0.01:tf, x0); % fixed step solver
%
% judge the intersection by extracting all temporal info
v1 = x(1:length(t)-1,index) - c;
v2 = x(2:length(t),index) - c;
ind = find(v1.*v2<0);
xtemp = x;
xtemp(:,index) = [];
%
if ~isempty(ind)
    ind = ind(1);  % only accept the first relocation
    if (all(xtemp(ind,:)'<ub) && all(xtemp(ind,:)'>lb)) && ...
            (all(xtemp(ind+1,:)'<ub) && all(xtemp(ind+1,:)'>lb))
        % ensure the image point is within the region of interest
        tj = t(ind);
        tj1 = t(ind+1);
        xj = x(ind,index);
        xj1 = x(ind+1,index);
        %
        % actual intersecting time
        tc = (tj1*(xj-c) - tj*(xj1-c))/(xj-xj1);
        %
        % interpolation of other coordinates
        x(:,index) = []; % delete intersection coordinate
        py = x(ind,:)' + (tc-tj)/(tj1-tj)*(x(ind+1,:)-x(ind,:))';
        status = 1;
    else
        x(:,index) = [];
        py = x(end,:)'; % out of the intersection plane
        status = 0;
    end
else
    x(:,index) = [];
    py = x(end,:)'; % out of the intersection plane
    status = 0;
end