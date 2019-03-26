function fszc = stepsize(f, sopt, h)
    %Select stepsize control to be used
    switch sopt
        case 1
            %Original step size control initial step=1
            fszc = @(x, fx, v, df, t)armijo(f, x, fx, v, df, t, h, 1E-4);
        case 2
            %Step size control move only one cell at the time
            fszc = @(x, fx, v, df, t)szc(f, x, fx, v, df, t, h, 1E-4);
        case 3
            fszc = @(x, fx, v, df, t)szc2(f, x, fx, v, df, t, h, 1E-4);
        case 4
            fszc = @(x, fx, v, df, t)szc3(f, x, fx, v, df, t, h, 1E-4);
        case 5
            fszc = @(x, fx, v, df, t)szc4(f, x, fx, v, df, t, h, 1E-4);
        case 6
            fszc = @(x, fx, v, df, t)szc5(f, x, fx, v, df, t, h, 1E-4);
        case 7
            fszc = @(x, fx, v, df, t)szc6(f, x, fx, v, df, t, h, 1E-4);
        case 8
            fszc = @(x, fx, v, df, t)szc7(f, x, fx, v, df, t, h, 1e-4, 0.9);
    end        
end