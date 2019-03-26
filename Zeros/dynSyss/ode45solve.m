function [t, x] = ode45solve(f, tspan, x0)
%
% Fixed step size RK-4 ode solver
dt = tspan(2)-tspan(1);
xold = x0;
x = zeros(length(tspan),length(x0));
%
for i = 1:length(tspan)
    ti = tspan(i);
    x(i,:) = xold';
    k1 = f(ti, xold);
    k2 = f(ti+dt/2, xold+dt/2*k1);
    k3 = f(ti+dt/2, xold+dt/2*k2);
    k4 = f(ti+dt, xold+dt*k3);
    xnew = xold + dt/6*(k1+2*k2+2*k3+k4);
    xold = xnew;
end
%
t = tspan;    