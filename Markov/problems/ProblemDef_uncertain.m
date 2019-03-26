function [lb, ub, lb_par, ub_par, dsys] = ProblemDef_uncertain(No)
% -------------------------------------------------------------------------
% Problem definition of the global analysis with parameter uncertainity,
% the definition requires both the searching information from state space
% and the parameter space.
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
        % Autonomous Duffing with uncertain parameter
        b = 1; k = 1;
        lb = [-1.5; -1.5];
        ub = [1.5;1.5];
        lb_par = 1.15; % assume uncertainity comes at damping
        ub_par = 1.25;
        tf = 0.1;
        dsys = @(x,par)Auto_Duffing_uncertain(x,b,k,tf,par);
end