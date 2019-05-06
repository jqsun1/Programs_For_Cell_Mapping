function [lb, ub, dsys, ineq, f, fmop] = ProblemDef(No, N, dopt, sopt, opt)
% -------------------------------------------------------------------------
% Problem definition of attractor finding mission. The definition includes
% all necessary cell mapping information and discrete point to point
% mapping function handle.
%
% Input argument:
%      No:    Problem number (a bunch of benchmark dynamics are listed)
%
% Output arguments:
%    lb, ub:  Lower and upper bound of region of interest
%      N:     Cell space partition
%     dsys:   Discrete dynamical system, either difference equation or
%             Poincare mapping in pointwise space
%     fmop:   Absolute value of nonlinear equations f(x) which will be used
%             for GCM-MOP searching at the first phase
% -------------------------------------------------------------------------
%
switch No
    case 1
        % Newton's iteration (2i2o test)
        lb = [-5, -5]';
        ub = [5, 5]';
        f = @(x)[4*x(1)*(x(1)^2+x(2)-11)+2*(x(1)+x(2)^2-7);
            2*(x(1)^2+x(2)-11)+4*x(2)*(x(1)+x(2)^2-7)];
        fmop = @(x)abs(f(x));
        df = @(x)[12*x(1)^2+4*x(2)-42, 4*x(1)+4*x(2);
            4*x(1)+4*x(2), 12*x(1)^2+4*x(1)-26];
        %
        % step size control
        h = (ub - lb)./N;
        fszc = stepsize(f, sopt, h);
        %
        % adaptive newton's search
        dsys = dynamic_systems(f, df, fszc, h, opt, dopt);
        ineq = [];
    case 2
        % Newton's iteration (cross test)
        lb = [-5, -5]';
        ub = [5, 5]';
        f = @(x)x(1)^2 - x(2)^2;
        fmop = @(x)abs(f(x));
        df = @(x)[2*x(1), -2*x(2)];
        %
        % step size control
        h = (ub - lb)./N;
        fszc = stepsize(f, sopt, h);
        %
        % adaptive newton's search
        dsys = dynamic_systems(f, df, fszc, h, opt, dopt);
        ineq = [];
    case 3
        % Newton's iteration (zig-zag saddle surface)
        lb = [-1, -1]';
        ub = [1, 1]';
        f = @(x)sin(0.5*x(1)^2-0.25*x(2)^2+3)*cos(2*x(1)+1-exp(x(2)));
        fmop = @(x)abs(f(x));
        df = @(x)[cos(0.5*x(1)^2-0.25*x(2)^2+3)*x(1)*cos(2*x(1)+1-exp(x(2)))-...
            2*sin(0.5*x(1)^2-0.25*x(2)^2+3)*sin(2*x(1)+1-exp(x(2))),...
            -0.5*x(2)*cos(0.5*x(1)^2-0.25*x(2)^2+3)*cos(2*x(1)+1-exp(x(2)))+...
            exp(x(2))*sin(0.5*x(1)^2-0.25*x(2)^2+3)*sin(2*x(1)+1-exp(x(2)))];
        %
        % step size control
        h = (ub - lb)./N;
        fszc = stepsize(f, sopt, h);
        %
        % adaptive newton's search
        dsys = dynamic_systems(f, df, fszc, h, opt, dopt);
        ineq = [];
    case 4
        % PD time-delay control system's stability region [kp, kd]
        lb = [-6; -5];
        ub = [5; 5];
        tau = pi/3; % time delay
        h = (ub - lb)./N;
        f = @(x)Test_PD_Stability(x, tau); % Maximum eigenvalue of CTA
        fmop = @(x)abs(f(x));
        df = @(x)cell_jacobian(f, x, h);
        %
        % step size control
        fszc = stepsize(f, sopt, h);
        %
        % adaptive newton's search
        dsys = dynamic_systems(f, df, fszc, h, opt, dopt);
        ineq = [];
    case 5
        lb = [-5; -5];
        ub = [5; 5];
        f = @(x)[x(1)^3-3*x(1)*x(2)^2-x(1)+1/sqrt(2);
            -x(2)^3 + 3*x(1)^2*x(2) - x(2)];
        fmop = @(x)abs(f(x));
        df = @(x)[3*x(1)^2-3*x(2)^2-1, 6*x(1)*x(2);
            6*x(1)*x(2), -3*x(2)^2+3*x(1)^2-1];
        h = (ub - lb)./N;
        fszc = stepsize(f, sopt, h);
        %
        % adaptive newton's search
        dsys = dynamic_systems(f, df, fszc, h, opt, dopt);
        ineq = [];
    case 6
        % high dim test problem
        n = length(N);
        lb = -0.3*ones(length(N),1);
        ub = 0.8*ones(length(N),1);
        f = @(x)More(x, n);
        fmop = @(x)abs(f(x));
        h = (ub - lb)./N;
%         df = @(x)cell_jacobian(f, x, h);
        df = @(x)Df_More(x);
        fszc = stepsize(f, sopt, h);
        %
        % adaptive newton's search
        dsys = dynamic_systems(f, df, fszc, h, opt, dopt);
        ineq = [];
    case 7
        lb = [-5;-3];
        ub = [5;3];
        f = @(x)[x(1)-2*x(2);x(1)*x(2)+x(1)-4*x(2)-4];
        fmop = @(x)abs(f(x));
        df = @(x)[1, -2;x(2)+1, x(1)-4];
        h = (ub - lb)./N;
        fszc = stepsize(f, sopt, h);
        dsys = dynamic_systems(f, df, fszc, h, opt, dopt);
        ineq = [];
    case 8
        % cluster case taken from oliver's paper
        lb = [-3;-3];
        ub = [3;3];
        h = (ub - lb)./N;
        f = @(x)fcluster(x);
        fmop = @(x)abs(f(x));
        df = @(x)cell_jacobian(f, x, h);
        fszc = stepsize(f, sopt, h);
        dsys = dynamic_systems(f, df, fszc, h, opt, dopt);
        ineq = [];
    case 9
        % high-dim polynominal
        n = length(N);
        lb = -15*ones(length(N),1);
        ub = 15*ones(length(N),1);
        h = (ub - lb)./N;
        lambda = 100;
        f = @(x)fpoly(x,lambda,n);
        fmop = @(x)abs(f(x));
        df = @(x)dfpoly(x,lambda,n);
        fszc = stepsize(f, sopt, h);
        dsys = dynamic_systems(f, df, fszc, h, opt, dopt);
        ineq = [];
    case 10
        % 2dof robotic PD control with time delay
        lb = [0; 0; -1.5; 0];
        ub = [10; 10; 3.5; 2.5];
        h = (ub - lb)./N;
        tau = pi/20;
        f = @(x)dof2CTApd(x, tau);
        fmop = @(x)abs(f(x));
        df = @(x)cell_jacobian(f, x, h);
        fszc = stepsize(f, sopt, h);
        dsys = dynamic_systems(f, df, fszc, h, opt, dopt);
        ineq = [];
    case 11
        % 2d workspace of a mechanism
        lb = [-2 -8 0 -pi/2 -pi/2 0 0]';
        ub = [8 8 pi/4 pi/2 pi/2 1 1]';
        h = (ub - lb)./N;
        f = @(x)mechanism(x);
        fmop = @(x)abs(f(x));
        df = @(x)Df_mechanism(f, x);
        fszc = stepsize(f, sopt, h);
        dsys = dynamic_systems(f, df, fszc, h, opt, dopt);
        ineq = [];
    case 12
        % chemical equilibrium application
        lb = [-0.5 40 -0.5 -0.5 -0.5]';
        ub = [0.5 50 0.5 0.5 0.5]';
        h = (ub - lb)./N;
        f = @(x)chemical(x);
        fmop = @(x)abs(f(x));
        df = @(x)Df_chemical(f,x);
        fszc = stepsize(f, sopt, h);
        dsys = dynamic_systems(f, df, fszc, h, opt, dopt);
        ineq = [];
    case 13
        % bio-example with time-delay (find all fixed points)
        S = 1; k1 = 0.05; Kd = 1; p = 5; kdx = 0.05;
        ksy = 1; kdy = 0.05; k2 = 1; ET = 1; Km = 0.1; KI = 2;
        lb = [0;0];
        ub = [1;2.5];
        h = (ub - lb)./N;
        f = @(x)bioRNA(x,S,k1,Kd,p,kdx,ksy,kdy,k2,ET,Km,KI);
        fmop = @(x)abs(f(x));
        df = @(x)Df_bioRNA(x,S,k1,Kd,p,kdx,ksy,kdy,k2,ET,Km,KI);
        fszc = stepsize(f, sopt, h);
        dsys = dynamic_systems(f, df, fszc, h, opt, dopt);
        ineq = [];
    case 14
        % 2dof neural network stable solutions
        mu1 = 0.2; mu2 = -1.3; w = 5; eps = 0.01; v1 = 1.45; v2 = -0.95; rho = -0.65;
        w12 = 0.35; w21 = 0.68;
        f = @(x)neural(x,mu1,mu2,w,rho,w12,w21,v1,v2,eps);
        lb = [-5;-5];
        ub = [5;5];
        h = (ub - lb)./N;
        df = @(x)Df_neural(f,x);
        fszc = stepsize(f, sopt, h);
        dsys = dynamic_systems(f, df, fszc, h, opt, dopt);
        ineq = [];
    case 15
        % oliver's map
        f = @(x)oliver(x);
        lb = [-2;-2;-2];
        ub = [2;2;2];
        h = (ub - lb)./N;
        df = @(x)Df_oliver(x);
        fszc = stepsize(f, sopt, h);
        dsys = dynamic_systems(f, df, fszc, h, opt, dopt);
        ineq = [];
    case 16
        % TD mathieu PD stability boundary
        eps = 1; delta = 4; tau = pi/4; k = 20;
        f = @(x)SD_mathieu(x,delta,eps,tau,k)-1;
        lb = [-8;-5];
        ub = [8;1];
        h = (ub - lb)./N;
        df = @(x)Df_SD_mathieu(f,x);
        fszc = stepsize(f, sopt, h);
        dsys = dynamic_systems(f, df, fszc, h, opt, dopt);
        ineq = [];
    case 17
        % Beale's function graident
         f = @(x)[(1.5-x(1)+x(1)*x(2))*(-1+x(2))+(2.25-x(1)+...
             x(1)*x(2)*x(2))*(-1+x(2)*x(2))+(2.625-x(1)+x(1)*x(2)^3)*(-1+x(2)^3);...
             (1.5-x(1)+x(1)*x(2))*x(1)+(2.25-x(1)+x(1)*x(2)*x(2))...
             *2*x(1)*x(2)+(2.625-x(1)+x(1)*x(2)^3)*3*x(1)*x(2)*x(2)];
        lb = [1;0];
        ub = [4;1];
        h = (ub - lb)./N;
        df = @(x)Df_SD_mathieu(f,x);
        fszc = stepsize(f, sopt, h);
        dsys = dynamic_systems(f, df, fszc, h, opt, dopt);
        ineq = [];
    case 18
        % camel function
        f = @(x)[4*x(1)-4.2*x(1)^3+x(1)^5+x(2);
            x(1)+2*x(2)];
        lb = [-2;-2];
        ub = [pi;pi];
        h = (ub - lb)./N;
        df = @(x)Df_SD_mathieu(f,x);
        fszc = stepsize(f, sopt, h);
        dsys = dynamic_systems(f, df, fszc, h, opt, dopt);
        ineq = [];
    case 19
        % easom function
        f = @(x)[-cos(x(2))*exp(-(x(1)-pi)^2-(x(2)-pi)^2)*(-sin(x(1))-2*(x(1)-pi)*cos(x(1)));
            -cos(x(1))*exp(-(x(1)-pi)^2-(x(2)-pi)^2)*(-sin(x(2))-2*(x(2)-pi)*cos(x(2)))];
        lb = [-4;-4];
        ub = [4;4];
        h = (ub - lb)./N;
        df = @(x)Df_SD_mathieu(f,x);
        fszc = stepsize(f, sopt, h);
        dsys = dynamic_systems(f, df, fszc, h, opt, dopt);
        ineq = [];
    case 20
        % fractional character equation
        f = @(x)[x(1)^sqrt(7)*cos(sqrt(7)*x(2))+3*x(1)^sqrt(5)*cos(sqrt(5)*x(2))+...
            5*x(1)^0.8*cos(0.8*x(2))+4;...
            x(1)^sqrt(7)*sin(sqrt(7)*x(2))+3*x(1)^sqrt(5)*sin(sqrt(5)*x(2))+...
            5*x(1)^0.8*sin(0.8*x(2))];
        lb = [0;0];
        ub = [30;2*pi];
        h = (ub - lb)./N;
        df = @(x)[sqrt(7)*x(1)^(sqrt(7)-1)*cos(sqrt(7)*x(2))+3*sqrt(5)*x(1)^(sqrt(5)-1)*cos(sqrt(5)*x(2))+4*x(1)^-0.2*cos(0.8*x(2)),...
            -x(1)^sqrt(7)*sqrt(7)*sin(x(2))-3*x(1)^sqrt(5)*sqrt(5)*sin(sqrt(5)*x(2))-4*x(1)^0.8*sin(0.8*x(2));...
            sqrt(7)*x(1)^(sqrt(7)-1)*sin(sqrt(7)*x(2))+3*sqrt(5)*x(1)^(sqrt(5)-1)*sin(sqrt(5)*x(2))+4*x(1)^-0.2*sin(0.8*x(2)),...
            x(1)^sqrt(7)*sqrt(7)*cos(x(2))+3*x(1)^sqrt(5)*sqrt(5)*cos(sqrt(5)*x(2))+4*x(1)^0.8*cos(0.8*x(2))];
        fszc = stepsize(f, sopt, h);
        dsys = dynamic_systems(f, df, fszc, h, opt, dopt);
        ineq = [];        
end