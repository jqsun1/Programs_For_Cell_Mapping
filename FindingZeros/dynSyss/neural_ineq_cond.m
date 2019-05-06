function indicator = neural_ineq_cond(mu1,mu2,w,rho,eps,h)
% conditions for the 2dof neural to have chaos, all conditions must be less
% than zero in order to satisify, rho must be (-1,0)
v = zeros(10,1);
v(1) = -(1-mu1-mu2)/w*eps;
v(2) = abs((1-mu1-mu2)/w) - abs(1/2/eps);
v(3) = (1-mu1-mu2)/w*eps + abs(h/w) - rho - 1;
v(4) = (1-mu1-mu2)/w*eps + abs(h/w) + rho;
v(5) = -eps;
v(6) = 1 - mu2 - w/2/eps;
v(7) = mu2 + 1;
v(8) = (-eps - w*(1+rho) - h)/mu2 + h - w*(1+rho) - mu2*eps;
v(9) = -mu2*eps + w*rho + h - (eps-w*rho+h)/mu2;
v(10) = 1 + abs(mu1) - mu2 - w/2/eps;
indicator = find(v<0);