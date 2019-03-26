function [t,x] = SD_mathieu_CTA(delta,eps,tau,tf,K,k)
kp = K(1);
kd = K(2);
A = @(t)[0,1;-(delta+2*eps*cos(2*t)),0];
Atau = @(t)[0,0;-kp,-kd];
n = length(A(0));
dt = pi/k; % discrete of the period
N = ceil(tau/dt); % discrete of time delay

% construct CTA matrix for time varying system
dtau = tau/N;
A_hat_subl = [1/dtau*eye(N*n), zeros(N*n,n)];
A_hat_subr = [zeros(N*n,n), 1/dtau*eye(N*n)];
A_hat_sup = @(t)[A(t), zeros(n,(N-1)*n), Atau(t)];
A_hat = @(t)[A_hat_sup(t);A_hat_subl-A_hat_subr];
x0 = zeros(n*(N+1),1); x0(1:n:end) = 0.1;

[t,x] = ode45(@(t,x)odemodel(t,x,A_hat),[0 tf],x0);

function dxdt = odemodel(t,x,A_hat)
dxdt = A_hat(t)*x;