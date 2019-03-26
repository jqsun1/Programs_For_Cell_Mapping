function MaxEig = Test_PD_Stability(v, tau)
% -------------------------------------------------------------------------
% Stability check via CTA approach, tracking equilibrum point at [r 0 0] 
% is examined.
% -------------------------------------------------------------------------
kp = v(1);
kd = v(2);
%
k = 4; c = 0.2;
N = 10;
A = [0 1;-k -c];
Atau = [0 0;-kp -kd];
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