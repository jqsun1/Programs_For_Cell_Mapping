function t = szc5(f, x, fx, v, df, t, h, c1)
    t = max(t);
    while all(f(x + t * v) >= fx)
        t = t/2;
        if t <= min(h)/2
            break;
        end
    end
end