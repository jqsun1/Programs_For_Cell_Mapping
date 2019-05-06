clc
clear
close all
PF = load('PF.dat');
PS = load('PS.dat');
figure(1)
scatter(PF(:,1),PF(:,2))
figure(2)
scatter(PS(:,1),PS(:,2))