function xnew = Newton(x, f, df, fszc, h, opt)
% classical newton's iteration
fx = f(x);
dx = df(x);
v = -pinv(dx)*fx;
% v = -dx\fx; % least square approximation
if opt
    v = v/norm(v);
end
%
t = fszc(x, fx, v, df, max(h));
%
xnew = x + t*v;