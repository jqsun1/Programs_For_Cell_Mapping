function x = bisect(f,lb,ub,tol,maxiter)
% -------------------------------------------------------------------------
% Bisection method for scalar nonlinear algebraic equation.
%
% By: Free Xiong; 2014-10-21
% -------------------------------------------------------------------------
if(nargin==3)
    tol = 1e-6;
    maxiter = 1000;
elseif(nargin==4)
    maxiter = 20000;
end

for i = 1:maxiter
    md = (ub+lb)/2;
    fmd = f(md);
    if(sign(f(lb))~=sign(fmd))
        ub = md;
    elseif(sign(f(ub))~=sign(fmd))
        lb = md;
    end
    if(ub-lb<tol)
        x = (ub+lb)/2;
        break;
    end
    x = (ub+lb)/2;
end