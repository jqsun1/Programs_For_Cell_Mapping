function dsys = dynamic_systems(f, df, fszc, h, opt, dopt)
% choose different local searching schemes
%
switch dopt
    case 1
        % classical newton iteration
        dsys = @(x)Newton(x, f, df, fszc, h, opt);
    case 2
        % modidifed newton iteration
        dsys = @(x)Newton_modified(x, f, df,  fszc, h, opt);
    case 3
        % broyden method without jacobian
        step = 15;
        dsys = @(x)Broyden(x, f, df, fszc, h, opt, step);
    case 4
        % trust region algorithm
        step = 1;
        dsys = @(x)Trust_region(x, f, df, fszc, h, opt, step);
end