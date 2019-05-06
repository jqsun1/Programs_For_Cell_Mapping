function res = neural(x,mu1,mu2,w,rho,w12,w21,v1,v2,eps)
res = zeros(2,1);
g1 = g(x(1),eps);
g2 = g(x(2),eps);
res(1) = (mu1+mu2-1)*x(1) + w*(g1+rho) + w12*g2 + v1;
res(2) = (mu1+mu2-1)*x(2) + w*(g2+rho) + w21*g1 + v2;
function g = g(x,eps)
% saturation function
g = 1/4*(2+abs(x/eps+1)-abs(x/eps-1));