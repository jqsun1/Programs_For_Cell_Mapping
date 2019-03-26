function tkp1 = multi_impact(tk,vk,u0,omega,omega1,u1)
% -------------------------------------------------------------------------
% Multiple impact condition of Eq (18) in Ref [1]
% Ref [1]: On Set-oriented Numerical Methods for 
%          Global Analysis of Non-smooth Mechanical Systems
% 
% If there is a t(k+1) solved, then the value is returned. Otherwise, -1 is
% returned as the indicator of zero multiple impact during one period.
%
% By: Free Xiong; 2015-01-27
% -------------------------------------------------------------------------
f = @(tkp1)condition(tkp1,tk,vk,u0,omega,omega1,u1);
dt = 2*pi/omega/100; % incremental interval
t0 = tk; % initial moment
tkp1 = tk+dt;
res_curr = f(tk);
res_next = f(tkp1);

while res_curr*res_next >= 0
    res_curr = res_next;
    tk = tkp1;
    tkp1 = tkp1+dt;
    res_next = f(tkp1);
    
    if res_curr*res_next < 0
        % tkp1 found
        tkp1 = bisect(f,tk,tkp1);
        break;
    elseif tkp1-t0 >= 2*pi/min([omega,omega1])
        % end of period reached, tkp1 not found
        tkp1 = -1;
        break;
    end
end

function res = condition(tkp1,tk,vk,u0,omega,omega1,u1)
res = u0/vk*(cos(omega*tkp1)-cos(omega*tk))-(tkp1-tk)+u1/vk*(sin(omega1*tkp1)-sin(omega1*tk));