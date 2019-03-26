function y = fpoly(x,lambda,n)
y = repmat(-sum(x),n,1)+lambda*x+x.^2-x.^3;