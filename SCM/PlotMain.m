clear all
close all
clc

% load data\duffing_nonauto_N_151x151.mat
% load data\vandepol_N_121x121.mat
% load data\henon_N_151x151.mat
% load data\duffing_noExcit_N_151x151.mat
% load data\duffing_auto_N_121x121.mat
% load data\gingerbreadman_N_151x151.mat
% load data\ikeda_N_151x151.mat
% load data\sun_N_151x151.mat
% load data\sun2_N_115x151.mat

xc = PlotCells(1:prod(N),lb,ub,N,'red','EdgeOff');
close

hold all
% plot attractors and different attraction domains in various colors
gru = setdiff(unique(gr),1);
c_att = linspace(0.4,1,length(gru));
c_dom = linspace(0.4,1,length(gru));

stc = 15; % the threshold step number for plot

% each ground is plotted with a slightly different color
for i = 1:length(gru)
    plot(xc(gr==gru(i) & st==0,1),xc(gr==gru(i) & st==0,2),'.',...
        'MarkerFaceColor',[c_att(i),0,0],...
        'MarkerEdgeColor',[c_att(i),0,0],...
        'MarkerSize',10);
    plot(xc(gr==gru(i) & st>stc,1),xc(gr==gru(i) & st>stc,2),'.',...
        'MarkerFaceColor',[0,0,c_dom(i)],...
        'MarkerEdgeColor',[0,0,c_dom(i)],...
        'MarkerSize',5);
end
axis([lb(1) ub(1) lb(2) ub(2)])
box on
% axis equal
xlabel('$x_1$')
ylabel('$x_2$')

% print -depsc duffing_nonauto_N_151x151.eps
% print -depsc vandepol_N_121x121.eps
% print -depsc henon_N_151x151.eps
% print -depsc duffing_noExcit_N_151x151.eps
% print -depsc duffing_auto_N_121x121.eps
% print -depsc gingerbreadman_N_151x151.eps
% print -depsc ikeda_N_151x151.eps
% print -depsc sun_N_151x151.eps
% print -depsc sun2_N_171x171.eps