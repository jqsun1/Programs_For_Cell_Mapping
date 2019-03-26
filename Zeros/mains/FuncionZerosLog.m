clear all
close all
clc
%
% log file that reports all necessary information of each example
load FunctionZeros_saddle.mat
disp('Saddle function example')
time
whos scm_att
Nm
iter
%
% log file that reports all necessary information of each example
load FunctionZeros_2I2O.mat
disp('2I2O function example')
time
whos scm_att
Nm
iter
%
load FunctionZeros_high_dim.mat
disp('4I4O function example')
time
whos scm_att
Nm
iter
%
load FunctionZeros_TD_PD.mat
disp('Time delay example')
time
whos scm_att
Nm
iter