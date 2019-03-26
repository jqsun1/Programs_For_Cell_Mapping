function t = szc(f, x, fx, v, df, t, h, c1)
    %compute default step size min(h)/abs(v)
    dx = df(x);
    t = min((t+10*eps) ./ (abs(v)));
    %Check for armijo condition and backtracking
    while all(f(x + t * v) > fx + c1 * t * dx * v)
        t = t/2;
        if t < min(h)/2
            break;
        end
    end
end