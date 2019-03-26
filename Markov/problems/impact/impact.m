function xkp1 = impact(xk,u0,omega,a1,a2,h)
% -------------------------------------------------------------------------
% Impact dynamics defiend in Eq (15) and (18) in Ref [1]. Note this
% dynamical system is defined in implicit algebrical form that requires
% careful numerical treatment.
%
% [1]  On Set-oriented Numerical Methods for
%      Global Analysis of Non-smooth Mechanical Systems
%
% By: Free Xiong; 2015-01-27
% -------------------------------------------------------------------------
% xk = [tk, vk];
tk = xk(1);
vk = xk(2);

% configuration for newton's method
tol = 1e-5;
funtol = 1e-4;
maxiter = 10000;

% mulitple impact condition
tkp1 = multi_impact(tk,vk,u0,omega);
if tkp1~=-1
    vkp1 = a1*a2*vk - (1+a1)*u0*omega*sin(omega*tkp1);
    xkp1 = [tkp1;vkp1];
else
    % solve Eq (15) numerically
    options = optimoptions('fsolve','Jacobian','on');
    f = @(xkp1)impact_dyn(xkp1,xk,u0,omega,a1,a2,h);
    xkp1 = fsolve(f,xk,options);
    
%     % newton's method is directly written here
%     x = xk;
%     dx = zeros(2,1);
%     for i = 1:maxiter
%         fx = f(x);
%         if norm(fx)<funtol
%             xkp1 = x;
%             break;
%         end
%         dx(2) = -fx(1)/(1-u0*omega/vk/a2*sin(omega*x(1)));
%         dx(1) = -fx(2)-dx(2)*(1+a1)*u0*omega^2*cos(omega*x(1));
%         if norm(dx)<tol
%             xkp1 = x;
%             break;
%         end
%         x = x+dx;
%     end
%     if i==maxiter
%         xkp1 = x;
%         warning('Maximum iteration number reached solving impact dynamics...');
%     end
end

% all timings are considered in one period
T = 2*pi/omega;
if xkp1(1)>T
    n = floor(xkp1(1)/T);
    xkp1(1) = xkp1(1)-n*T;
end