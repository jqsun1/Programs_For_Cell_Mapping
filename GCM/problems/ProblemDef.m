function [lb, ub, dsys] = ProblemDef(No)
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
%
switch No
    case 1
        % Duffing oscillator with external exicitation using time interval
        % sampling method to build the scm
        lb = [1.6; -3.5];
        ub = [3.6; 5.5];
        k = 0.25; B = 8.5; alpha = 0.05; tf = 2*pi;
        dsys = @(y)DuffingGraph(y, k, alpha, B, tf);
    case 2
        lb = [-8;-3;2];
        ub = [8;3;8];
        d = -6.02; r = 6.74; tf = 1;
        dsys = @(y)plasmaGraph(y,r,d,tf);
    case 3
        lb = [-1;-1];
        ub = [1;1];
        mu = 3.6; alpha = 0.025;
        dsys = @(x)logistic2d(x,mu,alpha);
    case 4
        % singular oscillator
        lb = [-1;-1];
        ub = [1;1];
        alpha = 1.05; sigma = 1.1; p = 0.012;
        q = 0.05; omega = 0.95;
        tf = pi/omega;
        dsys = @(x)AFM(x,alpha,sigma,p,q,omega,tf);
    case 5
        % Hsu's map
        lb = [-pi;-7];
        ub = [pi;7];
        alpha = 5.7; mu = 0.1*pi;
        c1 = (1-exp(-2*mu))/2/mu;
        d1 = exp(-2*mu);
        dsys = @(x)[x(1)-c1*alpha*sin(x(1))+c1*x(2);
            -d1*alpha*sin(x(1))+d1*x(2)];
    case 6
        % Autonomous Duffing system
        c = 1.15;
        lb = [-1.5;-1.5];
        ub = [1.5;1.5];
        tf = 0.1;
        dsys = @(x)Auto_Duffing(x,c,tf);
    case 7
        % Random duffing-van-de-pol system
        mu = 0.13;
        D1 = 0.013;
        D2 = 0.029;
        tf = pi;
        lb = [-3;-5];
        ub = [3;5];
        dsys = @(x)duffin_vanDePol_rand(x,mu,D1,D2,tf);
    case 8
        % Random van-de-pol system
        ksi = 0.05;
        eps = 1;
        Dw = 0.01;
        tf = 0.5;
        lb = [-4.5;-4.5];
        ub = [4.5;4.5];
        dsys = @(x)vanDePol_rand(x,ksi,eps,Dw,tf);
    case 9
        % 1d random system
        alpha = -1;
        beta = 2;
        Dw = 0.1;
        tf = 0.15;
        lb = -2;
        ub = 2;
        dsys = @(x)rand_1d(x,alpha,beta,Dw,tf);
    case 10
        % 2d impact study
        a1 = 0.8; a2 = 0.95; u0 = 0.125; omega = pi/2; h = 0.15; omega1 = pi/2; u1 = -0.05;
        % x = [tk;vk]
        lb = [0;-0.1];
        ub = [2*pi/min([omega,omega1]);0.6];
        dsys = @(xk)impact2(xk,u0,omega,a1,a2,h,omega1,u1);
end