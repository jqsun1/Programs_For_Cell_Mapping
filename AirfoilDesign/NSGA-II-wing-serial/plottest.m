clc
clear
close all
% for i = 0 : 30 : 100
% 	PF = load(strcat('GAPF',num2str(i),'.dat'));
% 	scatter(PF(:,1),PF(:,2),'MarkerFaceColor',rand(1,3));
% 	hold on
% end
ref = load('NSGA2_airfoil_R0.dat');
scatter3(ref(:,1),ref(:,2),ref(:,3))
%scatter(ref(:,1),ref(:,2))