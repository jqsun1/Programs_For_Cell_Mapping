function MaxEig = Test_PD_Stability_tracking(v, tau)
kp = v(1); kd = v(2); ki = v(3);
%
A = [0 1 0;0 0 0;-1 0 0];
Atau = [0 0 0;kp -kd -ki;0 0 0];
N = 10;
%
% construct CTA matrix and vector
n = length(A);
dtau = tau/N;
A_hat_subl = [1/dtau*eye(N*n), zeros(N*n,n)];
A_hat_subr = [zeros(N*n,n), 1/dtau*eye(N*n)];
A_hat_sup = [A, zeros(n,(N-1)*n), Atau];
A_hat = [A_hat_sup;A_hat_subl-A_hat_subr];
%
J = A_hat;
MaxEig = max(real(eig(J)));