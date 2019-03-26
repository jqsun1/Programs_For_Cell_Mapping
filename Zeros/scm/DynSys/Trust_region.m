function x_new = Trust_region(x, f, df, fszc, h, opt, step)
% Trust region algorithm with several step iteration
%
r_max = 2*max(h); % maximum trust region radius
zeta = 0.15;
r_old = 0.75*r_max;
x_old = x;
%
for i = 1:step
    %
    % calculate direction in the current trust region
    dx = df(x_old);
    fx = f(x_old);
    res = @(x)1/2*norm(f(x));
    m = @(v)res(x_old) + v'*dx*fx + 1/2*v'*(dx'*dx)*v; % model function
    v0 = -pinv(dx)*fx;
    v0 = r_old*v0/norm(v0);
    %  
    v = fmincon(m,v0,[],[],[],[],-inf,inf,@(x)constraint(x,r_old));
    rho = ( norm(fx)^2 - norm(f(x_old+v))^2 )/...
        ( norm(fx)^2 - norm(fx+dx*v)^2 );
    %
    % calculate next trust region radius
    if rho < 1/4
        r_new = 1/4*r_old;
    else
        if rho > 3/4
            r_new = min([2*r_old,r_max]);
        else
            r_new = r_old;
        end
    end
    %
    % iterate to the next step with new center
    if rho > zeta
        x_new = x_old + v;
    else
        x_new = x_old;
        break
    end
    %
    % update center if not break
    x_old = x_new;
    r_old = r_new;
end