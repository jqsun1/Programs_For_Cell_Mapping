function x_new = Henon(x_old, a, b)
x_new = zeros(2,1);
x_new(1) = 1-a*x_old(1)^2+x_old(2);
x_new(2) = b*x_old(1);