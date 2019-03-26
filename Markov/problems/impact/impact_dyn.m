function [res, J] = impact_dyn(xkp1,xk,u0,omega,a1,a2,h)
% -------------------------------------------------------------------------
% Impact dynamics in Eq (15).
%
% By: Free Xiong; 2015-01-27
% -------------------------------------------------------------------------
vk = xk(2); tk = xk(1);
vkp1 = xkp1(2); tkp1 = xkp1(1);
res = zeros(2,1);

res(1) = tkp1-tk-1/vk*(h*(1+1/a2)-u0*cos(omega*tk)-u0/a2*cos(omega*tkp1));
res(2) = vkp1-a1*a2*vk+(1+a1)*u0*omega*sin(omega*tkp1);

J = [1-u0*omega/vk/a2*sin(omega*tkp1),0;(1+a1)*u0*omega^2*cos(omega*tkp1),1];