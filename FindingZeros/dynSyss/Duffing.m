function dxdt = Duffing(t,x,k,alpha,B)
dxdt = zeros(2,1);
dxdt(1) = x(2);
dxdt(2) = -k*x(2) - alpha*x(1) - x(1)^3 + B*cos(t);