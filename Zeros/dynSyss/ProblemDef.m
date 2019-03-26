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
    case 3
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
    case 4
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
    case 5
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
    case 6
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
end