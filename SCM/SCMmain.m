clear all
close all
clc

% -------------------------------------------------------------------------
% SCM for nonliear dynamical systems without subdivision, global searching
% with transient information included. All examples presented in this code
% set are at 2d state space.
%
% By: Furui Xiong; 10/24/2015
% -------------------------------------------------------------------------

tic;
addpath(genpath('./')) % Add all the sub-folders of the current directory

% define problem, the short description of each problem is listed here:
% 1, duffing oscillator with external excitation
% 2, van-de-pol equation
% 3, henon map, https://en.wikipedia.org/wiki/H%C3%A9non_map
% 4, duffing oscillator without excitation
% 5, discrete duffing map, https://en.wikipedia.org/wiki/Duffing_map
% 6, discrete gingerbreadman map, https://en.wikipedia.org/wiki/Gingerbreadman_map
% 7, discrete ikeda map, https://en.wikipedia.org/wiki/Ikeda_map
% 8, discrete sun's map, from "Effects of small random uncertainties on non-linear systems studied by the generalized cell mapping method"
% 9, another example from the same reference above
prob_No = 10;
[lb, ub, N, dsys] = ProblemDef(prob_No);

% construct scm from ode/discrete point mapping
[C, S] = graph(dsys, ub, lb, N);

% unravel scm with sequential algorithm
[gr, pe, st] = search(C,N);

% % backward sorting of scm
% [ipg, igr] = backward_scm(S,C);

% plot attrator
att = S(gr~=1 & st==0);
whos att
time = toc
xc = PlotCells(att,lb,ub,N,'red','EdgeOn');

% save data\duffing_nonauto_N_151x151.mat
% save data\vandepol_N_121x121.mat
% save data\henon_N_151x151.mat
% save data\duffing_noExcit_N_151x151.mat
% save data\duffing_auto_N_121x121.mat
% save data\gingerbreadman_N_151x151.mat
% save data\ikeda_N_151x151.mat
% save data\sun_N_151x151.mat
% save data\sun2_N_115x151.mat