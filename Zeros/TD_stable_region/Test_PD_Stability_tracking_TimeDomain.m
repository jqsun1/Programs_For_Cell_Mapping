function [t,x] = Test_PD_Stability_tracking_TimeDomain(v,tau)
wn = 5; zeta = 0.1;
kp = v(1); kd = v(2);
eps = 0.1; r = 1.0;
%
A = [0 1;-wn^2 -2*zeta*wn^2];
Atau = [0 0;kp -kd];
N = 20;
%
% construct CTA matrix and vector
n = length(A);
dtau = tau/N;
A_hat_subl = [1/dtau*eye(N*n), zeros(N*n,n)];
A_hat_subr = [zeros(N*n,n), 1/dtau*eye(N*n)];
A_hat_sup = [A, zeros(n,(N-1)*n), Atau];
A_hat = [A_hat_sup;A_hat_subl-A_hat_subr];
%
tf = 50;
x0 = zeros(length(A_hat),1);
[t,x] = ode45(@(t,x)odefun(t,x,A_hat,eps,kp,r),[0 tf],x0);

function dxdt = odefun(t,x,A_hat,eps,kp,r)
F = zeros(length(A_hat),1); 
F(2) = -kp*r - eps*x(1)^3;
dxdt = A_hat*x + F;