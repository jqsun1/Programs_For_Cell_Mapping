clear all
close all
clc
%
% test of the backward searching program
wq = cputime;
load 'DuffingAttractorIllustration.mat'
%
[lb, ub, dyn_scm] = ProblemDef(3, N_new);
k = 0.25; B = 8.5; alpha = 0.02;
dyn = @(t,x)Duffing(t,x,k,alpha,B); % duffing example
%
tf = 2*pi;
[iter, V] = attraction_basin(S_fine, N_new, lb, ub, dyn_scm);
%
time = cputime - wq
%
PlotCells(V, lb, ub, N_new)
hold on
PlotCells(S_fine, lb, ub, N_new, 'green')
