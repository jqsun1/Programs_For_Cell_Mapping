function J = Df_bioRNA(x,S,k1,Kd,p,kdx,ksy,kdy,k2,ET,Km,KI)
J = zeros(2,2);
J(1,1) = -kdx;
J(1,2) = -k1*S*Kd^p*p*x(2)^(p-1)/(Kd^p+x(2)^p)^2;
J(2,1) = ksy;
J(2,2) = -kdy-k2*ET*(Km+x(2)+KI*x(2)^2-x(2)*(1+2*KI*x(2)))/(Km+x(2)+KI*x(2)^2)^2;