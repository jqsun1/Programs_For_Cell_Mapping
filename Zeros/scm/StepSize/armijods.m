%Input:
%   f:  Vector of functions
%   x:  Point in R^N
%   fx: f(x)
%   v:  Direction
%   a:  Vector of weights
%   t:  Initial step size
%   c1: Constant

function t = armijods(f, x, fx, v, a, t, c1)
    %Compute default step size
    t = min((t + 10*eps) ./ (abs(v)));    
    %Check armijo condition and backtracking
    while f(x + t * v)' * a > fx' * a - c1 * t * norm(a)^2
        t = t/2;
        if t < eps
            break;
        end
    end    
end