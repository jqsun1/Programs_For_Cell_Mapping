function y = mechanism(X)
% x = [u v ksi], where u is 2-by-1, z is 3-by-1 and ksi is 2-by-1. The
% problem has a 7i6o structure
x = X(1); y = X(2); v = X(3:5); ksi = X(6:7);

phi = [x - cos(pi/3*(sin(v(1))+sin(v(2))+sin(v(3))))-...
    2*cos(pi/3*(sin(v(1))+sin(v(2))))-4*cos(pi/3*sin(v(1)));...
    y - sin(pi/3*(sin(v(1))+sin(v(2))+sin(v(3))))-...
    2*sin(pi/3*(sin(v(1))+sin(v(2))))-4*sin(pi/3*sin(v(1)))];
y(1:2) = phi;

J = zeros(2,3);
J(1,1) = pi/3*(sin(pi/3*(sin(v(1))+sin(v(2))+sin(v(3))))*cos(v(1))+...
    2*sin(pi/3*(sin(v(1))+sin(v(2))))*cos(v(1)) + 4*sin(pi/3*sin(v(1)))*cos(v(1)));
J(1,2) = pi/3*(sin(pi/3*(sin(v(1))+sin(v(2))+sin(v(3))))*cos(v(2)) + ...
    2*sin(pi/3*(sin(v(1))+sin(v(2))))*cos(v(2)));
J(1,3) = pi/3*(sin(pi/3*(sin(v(1))+sin(v(2))+sin(v(3))))*cos(v(3)));
J(2,1) = pi/3*(-cos(pi/3*(sin(v(1))+sin(v(2))+sin(v(3))))*cos(v(1)) - ...
    2*cos(pi/3*(sin(v(1))+sin(v(2))))*cos(v(1)) - 4*cos(pi/3*sin(v(1)))*cos(v(1)));
J(2,2) = pi/3*(-cos(pi/3*(sin(v(1))+sin(v(2))+sin(v(3))))*cos(v(2)) - ...
    2*cos(pi/3*(sin(v(1))+sin(v(2))))*cos(v(2)));
J(2,3) = pi/3*(-cos(pi/3*(sin(v(1))+sin(v(2))+sin(v(3))))*cos(v(3)));
y(3:5) = J'*ksi;

y(6) = ksi(1)^2 + ksi(2)^2 - 1;

y = reshape(y,6,1);