function y = chemical(x)
% chemical equilibrum application in Crina Grosan's paper on MOP zeros
R = 10;
R5  =0.193;
R6 = 0.002597/sqrt(40);
R7 = 0.003448/sqrt(40);
R8 = 0.00001799/sqrt(40);
R9 = 0.0002155/sqrt(40);
R10 = 0.00003846/40;

y = [x(1)*x(2) + x(1) - 3*x(5);
     2*x(1)*x(2) + x(1) + x(2)*x(3)^2 + R8*x(2)-R*x(5)+...
     2*R10*x(2)^2 + R7*x(2)*x(3) + R9*x(2)*x(4);
     2*x(2)*x(3)^2 + 2*R5*x(3)^2 - 8*x(5) + R6*x(3) + R7*x(2)*x(3);
     R9*x(2)*x(4) + 2*x(4)^2 - 4*R*x(5);
     x(1)*(x(2)+1) + R10*x(2)^2 + x(2)*x(3)^2 + R8*x(2) + ...
     R5*x(3)^2 + x(4)^2 - 1 + R6*x(3) + R7*x(2)*x(3) + R9*x(2)*x(4)];