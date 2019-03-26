function res = bioRNA(x,S,k1,Kd,p,kdx,ksy,kdy,k2,ET,Km,KI)
res = zeros(2,1);
res(1) = k1*S*Kd^p/(Kd^p+x(end)^p) - kdx*x(1);
res(2) = ksy*x(1) - kdy*x(2) - k2*ET*x(2)/(Km+x(2)+KI*x(2)^2);