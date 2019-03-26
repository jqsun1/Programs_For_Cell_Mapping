function t = szc7(f, x, fx, v, df, t, h, c1, c2)
% Wolfe conditions
while ~all( f(x + t*v) <= fx + c1*t*df(x)*v & ...
            df(x + t*v)*v >= c2*df(x)*v )
    t = t/2;
    if t <= min(h)/2
        break;
    end
end