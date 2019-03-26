% plot gcm info of the duffing example
% by: Free Xiong; 11/09/2015
clear all
close all
clc

addpath(genpath('./')) % Add all the sub-folders of the current directory

load data/duffing_tran_N_100x100.mat
% load data/impact_N_189x189.mat

% plot invariant set cells
PlotCells(find(pg==1 & gr==1),lb,ub,N,'red','EdgeOff');
hold on
PlotCells(find(pg==1 & gr~=1),lb,ub,N,'blue','EdgeOff');

% plot attraction domain
PlotCells(find(Dm(:,2)==2),lb,ub,N,'green','EdgeOff');
xlabel('$x_1$')
ylabel('$x_2$')