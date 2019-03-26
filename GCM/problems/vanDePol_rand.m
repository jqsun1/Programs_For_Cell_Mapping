function py = vanDePol_rand(px,ksi,eps,Dw,tf)
f = @(t,x)odefun(t,x,ksi,eps,Dw);
[~,x] = ode45solve(f,0:tf/50:tf,px);
py = x(end,:)';

function dxdt = odefun(t,x,ksi,eps,Dw)
dxdt = zeros(2,1);
dxdt(1) = x(2);
dxdt(2) = -2*ksi*(eps*x(1)^2-1)*x(2)-x(1)+sqrt(Dw)*randn();