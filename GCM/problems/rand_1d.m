function py = rand_1d(px,alpha,beta,Dw,tf)
f = @(t,x)odefun(t,x,alpha,beta,Dw);
[~,x] = ode45solve(f,0:tf/50:tf,px);
py = x(end,:)';

function dxdt = odefun(t,x,alpha,beta,Dw)
dxdt = -alpha*x(1)-beta*x(1)^3+sqrt(Dw)*randn();