function t = szc4(f, x, fx, v, df, t, h, c1)
    t = max(t);
    dx = df(x);
    while all(f(x + t * v) > fx + c1 * t * dx * v)
        t = t/2;
        if t <= min(h)/2
            break;
        end
    end
    vn = floor(abs(v) + 1/2);
    t = diag(t ./ abs(v + (v==0)) .* vn);
end