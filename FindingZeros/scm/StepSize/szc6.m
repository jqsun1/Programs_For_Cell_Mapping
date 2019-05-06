function t = szc6(f, x, fx, v, df, t, h, c1)
    t = max(t);    
    while ~all( abs(fx) > abs(f(x+v*t)) )
        t = t/2;
        if t <= min(h)/2
            break;
        end
    end
end