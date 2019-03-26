clear all
close all
clc

% plot from txt files
fileID = fopen('ps.dat');
ps = textscan(fileID,'%f %f %f');
ps = cell2mat(ps);
fclose(fileID);
figure(1)
plot3(ps(:,1),ps(:,2),ps(:,3),'.')
% plot(ps(:,1),ps(:,2),'.')

fileID = fopen('pf.dat');
pf = textscan(fileID,'%f %f');
pf = cell2mat(pf);
fclose(fileID);
figure(2)
plot(pf(:,1),pf(:,2),'.')
% plot3(pf(:,1),pf(:,2),pf(:,3),'.')