function xnew = DynSys(xold, mu)
% point to point mapping with difference equation
%
xnew = zeros(2,1);
xnew(1) = (1-mu)*xold(2) + (2-2*mu+mu^2)*xold(1)^2;
xnew(2) = -(1-mu)*xold(1);