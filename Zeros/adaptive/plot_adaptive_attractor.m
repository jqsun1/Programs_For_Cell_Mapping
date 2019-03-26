clear all
close all
clc
%
load DuffingAttractor_adaptive_cell_size.mat
%
figure(1)
for i = 1:length(Q)
    S = Q(i).members;
    N = Q(i).division;
    PlotCells(S, lb, ub, N);
    hold on
end
%
% print -depsc DuffingAttractor_adaptive.eps