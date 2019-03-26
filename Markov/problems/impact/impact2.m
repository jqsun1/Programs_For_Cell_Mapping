function xkp1 = impact2(xk,u0,omega,a1,a2,h,omega1,u1)
% -------------------------------------------------------------------------
% Impact dynamics defiend in Eq (15) and (18) in Ref [1]. Note this
% dynamical system is defined in implicit algebrical form that requires
% careful numerical treatment.
%
% [1]  On Set-oriented Numerical Methods for
%      Global Analysis of Non-smooth Mechanical Systems
%
% By: Free Xiong; 2015-01-28
% -------------------------------------------------------------------------
% xk = [tk, vk];
tk = xk(1);
vk = xk(2);
t0 = tk;

% mulitple impact condition
tkp1 = multi_impact(tk,vk,u0,omega,omega1,u1);
if tkp1~=-1
    vkp1 = a1*a2*vk - (1+a1)*u0*omega*sin(omega*tkp1) ...
        + (1+a1)*u1*omega1*cos(omega1*tkp1);
    xkp1 = [tkp1;vkp1];
else
    % solve the first equation Eq (15) numerically using bisection
    f = @(tkp1)first_eqn(tkp1,tk,vk,h,a2,u0,omega,u1,omega1);
    dt = 2*pi/omega/100; % incremental interval
    tkp1 = tk+dt;
    res_curr = f(tk);
    res_next = f(tkp1);
    
    % bisection
    while res_curr*res_next >= 0
        res_curr = res_next;
        tk = tkp1;
        tkp1 = tkp1+dt;
        res_next = f(tkp1);
        
        if res_curr*res_next < 0
            % [tkp,tkp1] bracket found
            tkp1 = bisect(f,tk,tkp1);
            break;
        elseif tkp1-t0 >= 2*pi/min([omega,omega1])
            tkp1 = inf;
            break;
        end
    end
    
    % substitute to second equation of Eq(15)
    vkp1 = a1*a2*vk - (1+a1)*u0*omega*sin(omega*tkp1)...
        + (1+a1)*u1*omega1*cos(omega1*tkp1);
    xkp1 = [tkp1;vkp1];
end

% all timings are considered in one period
T = 2*pi/omega;
if xkp1(1)>T
    n = floor(xkp1(1)/T);
    xkp1(1) = xkp1(1)-n*T;
end

function res = first_eqn(tkp1,tk,vk,h,a2,u0,omega,u1,omega1)
% first equation of Eq(15)
res = tkp1-tk-1/vk*(h*(1+1/a2)-u0*cos(omega*tk)-u0/a2*cos(omega*tkp1)...
    -u1*sin(omega1*tk)-u1/a2*sin(omega1*tkp1));