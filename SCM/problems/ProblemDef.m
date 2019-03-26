function [lb, ub, N, dsys] = ProblemDef(No)
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
%             Poincare mapping in pointwise space, for some cases, dsys is
%             a function handle with respect to 'y' and 'tf', where the
%             latter means the integration time for local mapping building
% -------------------------------------------------------------------------

switch No
    case 1
        % Duffing oscillator with external exicitation using time interval
        lb = [1.6; -3.5];
        ub = [3.6; 5.5];
        N = [151; 151];
        k = 0.25; B = 8.5; alpha = 0.05; tf = 2*pi;
        dsys = @(y)DuffingGraph(y, k, alpha, B, tf);
    case 2
        % Van-der-pol oscillator
        lb = [-3; -3];
        ub = [3; 3];
        mu = 0.25;
        N = [121; 121];
        dsys = @(y)VanDePol(y,mu);
    case 3
        % Henon map
        a = 1.4; b = 0.3;
        dsys = @(x)Henon(x, a, b);
        lb = [-2;-2];
        ub = [2;2];
        N = [151; 151];
    case 4
        % Autonomous Duffing system
        c = 1.2;
        lb = [-1.5;-1.5];
        ub = [1.5;1.5];
        N = [151; 151];
        tf = 1;
        dsys = @(x)Auto_Duffing(x,c,tf);
    case 5
        % Duffing map
        lb = [-2;-2];
        ub = [2;2];
        a = 2.75; b = 0.2;
        N = [121; 121];
        dsys = @(x)[x(2);-b*x(1)+a*x(2)-x(2)^3];
    case 6
        % Gingerbreadman map
        lb = [-14;-14];
        ub = [14;14];
        N = [151; 151];
        dsys = @(x)[1-x(2)+abs(x(1));x(1)];
    case 7
        % Ikeda map
        lb = [-5;-5];
        ub = [5;5];
        N = [151; 151];
        u = 0.8;
        t = @(x)0.4-6/(1+x(1)^2+x(2)^2);
        dsys = @(x)[1+u*(x(1)*cos(t(x))-x(2)*sin(t(x)));...
            u*(x(1)*sin(t(x))+x(2)*cos(t(x)))];
    case 8
        % Sun's map
        lb = [-1.5;-2.5];
        ub = [1.5;1.5];
        N = [151; 151];
        dsys = @(x)[0.9*x(2)+1.81*x(1)^2;-0.9*x(1)];
    case 9
        % Sun's 2nd map
        lb = [-pi;-7];
        ub = [pi;7];
%         lb = [0;-7];
%         ub = [pi;0];
        N = [115;151];
%         N = [191;211];
        mu = 0.1*pi;
        alpha = 5.7;
        c1 = (1-exp(-2*mu))/2/mu;
        d1 = exp(-2*mu);
        dsys = @(x)[x(1)-c1*alpha*sin(x(1))+c1*x(2);...
            -d1*alpha*sin(x(1))+d1*x(2)];
    case 10
        % 9 order vandepol
        lb = [-3;-3];
        ub = [3;3];
        N = [51;51];
        eps = -0.2;
        alpha = [2.5, 4.545, 2.5, 0.4];
        T = 4*pi;
        dsys = @(x)vandepol9(x,eps,alpha,T);
end