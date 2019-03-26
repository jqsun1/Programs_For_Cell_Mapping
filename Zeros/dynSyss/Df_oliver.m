function df = Df_oliver(x)
df = zeros(1,3);
df(1) = 4*x(1)^3 - 3*x(2) - sin(x(1)*x(2))*x(2);
df(2) = -3*x(1) - sin(x(1)*x(2))*x(1);
df(3) = 4*sin(4*x(3));