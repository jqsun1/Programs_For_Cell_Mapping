clear all
close all
clc
% Test the PoincareMap code with one specific example (Lorenz system here)
%
y = [-2.5 3.5]'; % user defined point on the intersection
index = 2;
c = 6*sqrt(2);
sigma = 10; beta = 8/3; rho = 28;
f = @(t,x)Lorenz(t,x,sigma,rho,beta);
lb = [-10; -5];
ub = [10; 30];
%
[py, status] = PoincareMap(y, index, c, f, lb, ub)