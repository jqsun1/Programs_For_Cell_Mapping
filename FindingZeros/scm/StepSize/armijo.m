%Input:
%   f:  Vector of functions
%   x:  Point in R^N
%   fx: f(x)
%   v:  Direction
%   a:  Vector of weights
%   t:  Initial step size
%   c1: Constant

function t = armijo(f, x, fx, v, df, t, h, c1)
    %Check for armijo condition
    dx = df(x);
    while all(f(x + t * v) >= fx + c1 * t * dx * v)
        t = t/2;
        if t <= min(h)/2
            break;
        end
    end
end