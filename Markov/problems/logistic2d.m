function y = logistic2d(x,mu,alpha)
y = zeros(2,1);
y(1) = (1-alpha)*mu*x(1)*(1-x(1))+alpha*mu*x(2)*(1-x(2));
y(2) = (1-alpha)*mu*x(2)*(1-x(2))+alpha*mu*x(1)*(1-x(1));