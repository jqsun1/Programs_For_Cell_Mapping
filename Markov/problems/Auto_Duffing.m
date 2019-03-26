function px = Auto_Duffing(y,c,tf)
% [~, x] = ode45(@(t,x)odefun(t,x,c),[0 tf],y);
[~, x] = ode45solve(@(t,x)odefun(t,x,c), 0:tf/20:tf, y);
px = x(end,:)';
%
function dxdt = odefun(t,x,c)
dxdt = zeros(2,1);
dxdt(1) = x(2);
dxdt(2) = -c*x(2)+x(1)-x(1)^3;