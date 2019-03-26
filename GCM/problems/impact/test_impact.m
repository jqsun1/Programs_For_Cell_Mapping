clear all
close all
clc

% test of impact dynamics code
a1 = 0.8; a2 = 0.9; u0 = 0.125; omega = pi/2; h = 0.15; omega1 = pi; u1 = 0.11;
xk = [0.1;0.1];
xkp1 = impact2(xk,u0,omega,a1,a2,h,omega1,u1)