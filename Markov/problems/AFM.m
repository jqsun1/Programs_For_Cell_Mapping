function py = AFM(y,alpha,sigma,p,q,omega,tf)
[~, x] = ode45solve(@(t,x)AFM_model(t,x,alpha,sigma,p,q,omega), 0:tf/100:tf, y);
py = x(end,:)';