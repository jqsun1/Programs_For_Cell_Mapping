clear all
close all
clc
%
% load FunctionZeros_TD_PD.mat
% load FunctionZeros_TD_PD_simplified.mat
% load FunctionZeros_saddle.mat
% load FunctionZeros_2I2O.mat
% load DuffingAttractorIllustration.mat
% load FunctionZeros_high_dim.mat
% load FunctionZeros_2I2O_N_13x9.mat
% load FunctionZeros_fcluster_GCM_SCM.mat
% load FunctionZeros_2dof_TD_PD_N_5x5x5x5_iter_2.mat
% load FunctionZeros_oliver3D_N_10x10x10_iter_2.mat

% Data used in the paper:
load FunctionZeros_2I2O_N_13x13_iter_4.mat
% load FunctionZeros_saddle_N_10x10_iter_4.mat
% load FunctionZeros_TD_PD_N_15x15_iter_4.mat
% load FunctionZeros_TD_PD_simplified_N_15x15_iter_4.mat   % only figure 3
% load FunctionZeros_poly_N_10x10x10x10_iter_3.mat
% load FunctionZeros_poly_N_15x15x15x15_iter_3.mat
% load FunctionZeros_TD_Mathieu_PD_N_10x10_iter_3.mat

figure(1)
xc1 = PlotCells(S, lb, ub, N_old, 'red', 'EdgeOn'); % scm sorting cells (coarse)
% print -depsc FuctionZeros_TD_PD_cover_set.eps
% print -depsc FunctionZeros_saddle_cover_set.eps
% print -depsc FunctionZeros_2I2O_cover_set.eps
% print -depsc FunctionZeros_2I2O_N_13x13_cover_set.eps
% print -depsc FunctionZeros_2I2O_N_13x9_cover_set.eps
% print -depsc FunctionZeros_TD_Mathieu_cover_set.eps

% figure(2)
% PlotCells(S1, lb, ub, N1, 'red', 'EdgeOn'); % 1st refined cells with block info avaliable
% print -depsc FunctionZeros_2I2O_N_13x13_FirstRefine.eps
% print -depsc FunctionZeros_2I2O_N_13x9_FirstRefine.eps

figure(3)
xc = PlotCells(S_new, lb, ub, N_new, 'red', 'EdgeOff'); % Refined attractor with expansion (fine)
% hold on
% plot(xa(:,1),xa(:,2),'b.') % only 2i2o and high-dim poly have fsolve solutions
% print -depsc FunctionZeros_TD_PD_final.eps
% print -depsc FunctionZeros_TD_PD_simplified_final.eps
% print -depsc FunctionZeros_saddle_final.eps
% print -depsc FunctionZeros_2I2O_final.eps
% print -depsc FunctionZeros_2I2O_final_comparison.eps
% print -depsc FunctionZeros_high_dim_final.eps
% print -depsc FunctionZeros_2I2O_N_13x13_final_miss.eps
% print -depsc FunctionZeros_2I2O_N_13x13_final_complete.eps
% print -depsc FunctionZeros_oliver3D_N_10x10x10_final.eps % view(60,10)
% print -depsc FunctionZeros_TD_Mathieru_final.eps

% hold on
% [X,Y] = meshgrid(-1:0.01:1);
% Z = sin(0.5*X.^2-0.25*Y.^2+3).*cos(2*X+1-exp(Y));
% mesh(X,Y,Z)
% hold on
% hold on
% surf(X,Y,zeros(size(Z)),'EdgeColor','None')
% axis([lb(1) ub(1) lb(2) ub(2)])
% print -depsc FunctionZeros_saddle_final_comparison.eps

% figure
% subplot(1,2,1)
% xc1 = PlotCells(S, lb, ub, N_old, 'red', 'EdgeOn'); % scm sorting cells (coarse)
% axis square
% subplot(1,2,2)
% xc = PlotCells(S_new, lb, ub, N_new, 'red', 'EdgeOff'); % Refined attractor with expansion (fine)
% axis square
% print -depsc FunctionZeros_2I2O_cover_set.eps