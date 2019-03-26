function x_new = Broyden(x, f, df, fszc, h, opt, step)
% Broyden's method with approximation of jacobian matrix, we need several
% steps in the local searching algorithm to approach a decent estimation of
% local jacobian matrix.
B = eye(size(df(x)));
x_old = x;
counter = 0;
%
while true
    v = -pinv(B)*f(x_old);
    t = fszc(x_old, f(x_old), v, df, 1);
    x_new = x_old + t*v;
    %
    if counter == step
        break
    else
        s = x_new - x_old;
        y = f(x_new) - f(x_old);
        %
        % update the so-called jacobian matrix
        B = B + (y-B*s)*s'/(s'*s);
        counter = counter + 1;
        x_old = x_new;
    end
end