function px = Auto_Duffing_uncertain(y,b,k,tf,par)
% [~, x] = ode45(@(t,x)odefun(t,x,c),[0 tf],y);
[~, x] = ode45solve(@(t,x)odefun(t,x,b,k,par), 0:tf/20:tf, y);
px = x(end,:)';
%
function dxdt = odefun(t,x,b,k,par)
dxdt = zeros(2,1);
dxdt(1) = x(2);
dxdt(2) = -par*x(2)+k*x(1)-b*x(1)^3;