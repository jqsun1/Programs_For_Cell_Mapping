function [t,x] = SD_mathieu_ode(delta,eps,tau,tf,K)
kp = K(1);
kd = K(2);
A = @(t)[0,1;-(delta+2*eps*cos(2*t)),0];
Atau = @(t)[0,0;-kp,-kd];
h = @(t)ddehist(t); % constant history
f = @(t,y,Z)ddemodel(t,y,Z,A,Atau); % dde model
sol = dde23(f,tau,h,[0 tf]);
t = sol.x;
x = (sol.y)';

function s = ddehist(t)
s = [0.1;0];

function dydt = ddemodel(t,y,Z,A,Atau)
ylag = Z(:,1);
dydt = A(t)*y + Atau(t)*ylag;