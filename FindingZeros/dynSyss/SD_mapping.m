function phi = SD_mapping(A,Atau,T,tau,k)
% -------------------------------------------------------------------------
% Determine the one-period mapping of linear time varying system subject
% to a constant delay. The stability boundary will be determined by examine
% the maximum eigenvalue which lies in the unit circle. Semi-discretization
% method is applied here for time delay approximation.
%
% Input arguments:
%     A, Atau:    function handles represent the linear and delay parts
%       T:        unique period of the system, if either A and Atau is
%       tau:      time delay reflected in matrix Atau
%       k:        discretization of the period T
%
% Output argument:
%      phi:       periodic mapping determines the stability boundary
% -------------------------------------------------------------------------
dt = T/k; % discrete of the period
N = ceil(tau/dt); % discrete of time delay
t = 0:dt:T;
n = length(A(0)); % dim of the system
%
% revolution a period for the mapping
phi = eye((N+1)*n);
for i = 1:length(t)-1
    Ai = A((t(i)+t(i+1))/2);
    Adi = Atau((t(i)+t(i+1))/2);
    Pi = dt*expm(Ai*(dt- (t(i)+t(i+1))/2 ))*Adi;
    Qi = expm(Ai*dt);
    Hi = [[Qi,zeros(n,(N-1)*n),Pi];[eye(n*N),zeros(N*n,n)]];
    phi = Hi*phi;
end