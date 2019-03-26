function J = dfpoly(x,lambda,n)
J = diag(1+(1-lambda)+2*x-3*x.^2)-ones(n);