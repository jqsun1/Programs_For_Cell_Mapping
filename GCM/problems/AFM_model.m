function dxdt = AFM_model(t,x,alpha,sigma,p,q,omega)
dxdt = zeros(2,1);
d = 4/27;
dxdt(1) = x(2);
dxdt(2) = -x(1)+d/(alpha-x(1))^2-sigma^6*d/30/(alpha-x(1))^8+...
    q*sin(omega*t)+p*q*omega*cos(omega*t)-p*x(2);