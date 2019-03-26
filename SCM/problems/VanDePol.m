function py = VanDePol(y, mu)
[~, x] = ode45solve(@(t,x)odefun(t,x,mu),0:2*pi/50:2*pi,y);
% [~, x] = ode45(@(t,x)odefun(t,x,mu),[0 1.3],y);
py = x(end,:)';
%
function dxdt = odefun(t,x,mu)
dxdt = zeros(2,1);
dxdt(1) = x(2);
dxdt(2) = mu*(1-x(1)^2)*x(2) - x(1);