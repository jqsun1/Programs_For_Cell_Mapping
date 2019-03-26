function xp = plasmaGraph(xc,r,d,tf)
f = @(t,x)plasma(t,x,r,d);
[~,x] = ode45solve(f,0:tf/100:tf,xc);
xp = x(end,:)';