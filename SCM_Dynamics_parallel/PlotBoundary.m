clear all
close all
clc
% plot boundary from scm whole domain sweeping
% By: Free Xiong; 2015-03-06

% load potential2d_zerofinding.mat
% load potential3d_zerofinding.mat

% mcells = load('boundary_cells.dat');
figure
if(length(lb)==2)
    xb = PlotCells(mcells, lb, ub, N, 'blue', 'EdgeOff');  
    hold on
    x = linspace(lb(1),ub(1));
    y = linspace(lb(2),ub(2));
    z = zeros(length(y),length(x));
    for i = 1:length(x)
        for j = 1:length(y)
            z(j,i) = f([x(i);y(j)]);
        end
    end
    contour(x,y,z,25);
    xlabel('$x$')
    ylabel('$y$')
    hold on
    plot(point_sol(sd_ind==0,1),point_sol(sd_ind==0,2),'r*') % sink
    hold on
    plot(point_sol(sd_ind==1,1),point_sol(sd_ind==1,2),'ro','MarkerFaceColor','r') % saddle
    hold on
    plot(point_sol(sd_ind==2,1),point_sol(sd_ind==2,2),'rx') % source  
end


if(length(lb)>=3)
    hold on
    x = PlotCells(mcells, lb, ub, N, 'blue', 'EdgeOff');
    close
    figure
    scatter3(x(:,1),x(:,2),x(:,3),1,abs(x(:,1)),'square','filled')
    hold on
    plot3(point_sol(sd_ind==0,1),point_sol(sd_ind==0,2),point_sol(sd_ind==0,3),'ko',...
        'MarkerFaceColor','k','MarkerSize',10) % sink
    hold on
    plot3(point_sol(sd_ind==1,1),point_sol(sd_ind==1,2),point_sol(sd_ind==1,3),'ro',...
        'MarkerFaceColor','r','MarkerSize',10) % saddle
    hold on
    plot3(point_sol(sd_ind==2,1),point_sol(sd_ind==2,2),point_sol(sd_ind==2,3),'rx') % source
    axis([lb(1) ub(1) lb(2) ub(2) lb(3) ub(3)])
    xlabel('$x_2$')
    ylabel('$x_3$')
    zlabel('$y_3$')
    view(15,26)
    print -depsc potential3d_full.eps
    
    % plot only half to display the boundary more clearly
    figure
    scatter3(x(x(:,1)>0,1),x(x(:,1)>0,2),x(x(:,1)>0,3),1,x(x(:,1)>0,1),'square','filled')
    hold on
    plot3(point_sol(sd_ind==0&point_sol(:,1)>0,1),...
        point_sol(sd_ind==0&point_sol(:,1)>0,2),...
        point_sol(sd_ind==0&point_sol(:,1)>0,3),'ko',...
        'MarkerFaceColor','k','MarkerSize',10) % sink
    hold on
    plot3(point_sol(sd_ind==1&point_sol(:,1)>0,1),...
        point_sol(sd_ind==1&point_sol(:,1)>0,2),...
        point_sol(sd_ind==1&point_sol(:,1)>0,3),'ro',...
        'MarkerFaceColor','r','MarkerSize',10) % saddle
    hold on
    plot3(point_sol(sd_ind==2&point_sol(:,1)>0,1),...
        point_sol(sd_ind==2&point_sol(:,1)>0,2),...
        point_sol(sd_ind==2&point_sol(:,1)>0,3),'rx') % source
    axis([0 ub(1) lb(2) ub(2) lb(3) ub(3)])
    xlabel('$x_2$')
    ylabel('$x_3$')
    zlabel('$y_3$')
    view(57,6)
    print -depsc potential3d_half_1.eps
    
    % half data plot from another perspective
    figure
    scatter3(x(x(:,1)>0,1),x(x(:,1)>0,2),x(x(:,1)>0,3),1,x(x(:,1)>0,1),'square','filled')
    hold on
    plot3(point_sol(sd_ind==0&point_sol(:,1)>0,1),...
        point_sol(sd_ind==0&point_sol(:,1)>0,2),...
        point_sol(sd_ind==0&point_sol(:,1)>0,3),'ko',...
        'MarkerFaceColor','k','MarkerSize',10) % sink
    hold on
    plot3(point_sol(sd_ind==1&point_sol(:,1)>0,1),...
        point_sol(sd_ind==1&point_sol(:,1)>0,2),...
        point_sol(sd_ind==1&point_sol(:,1)>0,3),'ro',...
        'MarkerFaceColor','r','MarkerSize',10) % saddle
    hold on
    plot3(point_sol(sd_ind==2&point_sol(:,1)>0,1),...
        point_sol(sd_ind==2&point_sol(:,1)>0,2),...
        point_sol(sd_ind==2&point_sol(:,1)>0,3),'rx') % source
    axis([0 ub(1) lb(2) ub(2) lb(3) ub(3)])
    xlabel('$x_2$')
    ylabel('$x_3$')
    zlabel('$y_3$')
    view(115,10)
    print -depsc potential3d_half_2.eps
    
    % main frame of boundary surface along x axis
    figure
    eps = 0.02;
    xd = linspace(eps,ub(1),18);
    xp = [];
    for i = 1:length(xd)
        xx = x(x(:,1)<xd(i)+eps & x(:,1)>xd(i)-eps,:);
%         plot3(xx(:,1),xx(:,2),xx(:,3),'.')
        xp = [xp; xx];
        hold on
    end
    scatter3(xp(:,1),xp(:,2),xp(:,3),2,xp(:,1),'filled','square')
    axis([0.1 ub(1) lb(2) ub(2) lb(3) ub(3)])
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$z$')
    hold on
    plot3(point_sol(sd_ind==0&point_sol(:,1)>0,1),...
        point_sol(sd_ind==0&point_sol(:,1)>0,2),...
        point_sol(sd_ind==0&point_sol(:,1)>0,3),'ko',...
        'MarkerFaceColor','k','MarkerSize',10) % sink
    hold on
    plot3(point_sol(sd_ind==1&point_sol(:,1)>0,1),...
        point_sol(sd_ind==1&point_sol(:,1)>0,2),...
        point_sol(sd_ind==1&point_sol(:,1)>0,3),'ro',...
        'MarkerFaceColor','r','MarkerSize',10) % saddle
    hold on
    plot3(point_sol(sd_ind==2&point_sol(:,1)>0,1),...
        point_sol(sd_ind==2&point_sol(:,1)>0,2),...
        point_sol(sd_ind==2&point_sol(:,1)>0,3),'rx') % source    
    view(151,8)
    grid on
    
    figure
    plot3(x(1:end,1),x(1:end,2),x(1:end,3),'b.','MarkerSize',2)
%     scatter3(x(1:end,1),x(1:end,2),x(1:end,3),2,-x(:,2),'filled','square')
    axis([0.1 ub(1) lb(2) ub(2) -1 2])
    view(116,14)
    xlabel('$x_2$')
    ylabel('$x_3$')
    zlabel('$y_3$')
    grid on
    print -depsc potential3d_half_open.eps
end