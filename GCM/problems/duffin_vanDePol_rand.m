function py = duffin_vanDePol_rand(px,mu,D1,D2,tf)
f = @(t,x)odefun(t,x,mu,D1,D2);
[~,x] = ode45solve(f,0:tf/100:tf,px);
py = x(end,:)';

function dxdt = odefun(t,x,mu,D1,D2)
dxdt = zeros(2,1);
a1 = 1.51;
a2 = 2.85;
a3 = 1.693;
a4 = 0.312;
dxdt(1) = x(2);
dxdt(2) = (-mu+a1*x(1)^2-a2*x(1)^4+a3*x(1)^6-a4*x(1)^8)*x(2)-x(1)-2*x(1)^3+...
    sqrt(D1)*randn()+x(1)*sqrt(D2)*randn();