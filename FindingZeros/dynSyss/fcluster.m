function g = fcluster(x)
g = zeros(2,1);
i = 1:20;
j = 21:40;
d = [0.1-0.01+i/1000, 0.9-0.03+j/1000];
g(1) = sin(4*x(2))*prod(x(1)-d);
g(2) = sin(4*x(1))*prod(x(2)-d);