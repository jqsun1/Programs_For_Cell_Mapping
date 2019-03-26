clear all
close all
clc

% test parametrized color plot
% By; Furui Xiong; 10/27/2015

% plot each column with different color
x = linspace(0,1);
y = linspace(0,1);
hold all

swEPSfigure
c = linspace(0.4,1,length(x));
for i = 1:length(x)
    scatter(x(i)*ones(size(y)),y,[],[c(i) 0 0],'filled','square')
end
axis equal
set(gca,'xlabel',[],'ylabel',[],'visible','off')