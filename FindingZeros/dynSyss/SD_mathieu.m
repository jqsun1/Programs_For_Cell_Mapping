function lambda = SD_mathieu(K,delta,eps,tau,k)
kp = K(1);
kd = K(2);
A = @(t)[0,1;-(delta+2*eps*cos(2*t)),0];
Atau = @(t)[0,0;-kp,-kd];
T = pi; % period of Mathieu equation
phi = SD_mapping(A,Atau,T,tau,k);
lambda = max(abs(eig(phi)));