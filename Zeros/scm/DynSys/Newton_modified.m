function xnew = Newton_modified(x, f, df, fszc, h, opt)
% modified newton's iteration
fx = f(x);
dx = df(x);
v = -pinv(dx)*fx;
if opt
    v = v/norm(v);
end
%
cos = -v'*(dx'*fx)/norm(v)/norm(dx'*fx);
delta = sqrt(3)/2;
%
if cos >= delta
    t = fszc(x, fx, v, df, 1.0);
    xnew = x + t*v;
else
    tau = 10;
    while true
        v = -inv( (dx'*dx+tau*eye(length(x))) )*dx'*fx;
        cos = -v'*(dx'*fx)/norm(v)/norm(dx'*fx);
        if cos >= delta
            break
        else
            tau = 2*tau;
        end
    end
    t = fszc(x, fx, v, df, 1.0);
    xnew = x + t*v;
end