function y = More(x, n)
I = diag(1:n);
y = n - repmat(sum(cos(x)),n,1) + I*(1-cos(x)) - sin(x); 