function dxdt = plasma(t,x,r,d)
dxdt = zeros(3,1);
dxdt(1) = x(1)+2*x(1)*x(2)^2*sin(x(3));
dxdt(2) = -r*x(2)-x(1)^2*x(2)*sin(x(3));
dxdt(3) = -2*d+2*(x(2)^2-x(1)^2)+2*(2*x(2)^2-x(1)^2)*cos(x(3));