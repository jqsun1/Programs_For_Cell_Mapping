function xc = PlotCells(S, lb, ub, N, color, edge, i, j, k)
% Visualize certain cell set S
%
if nargin < 5
    color = 'red';
    edge = 'EdgeOff';
elseif nargin < 7
    i = 1;
    j = 2;
    k = 3;
end

xc = zeros(length(S),length(N));
h = (ub - lb)./N;

swEPSfigure;
% swFigSize;

%
for index = 1: length(S)
    cell = S(index);
    z = celltoz(cell, N);
    x = ztox(z, h, lb);
    
    if strcmp(edge,'EdgeOn')
        if length(N) == 2
            % 2d square plot
            x_sw = [x(1) - h(1)/2; x(2) - h(2)/2];
            x_se = [x(1) + h(1)/2; x(2) - h(2)/2];
            x_nw = [x(1) - h(1)/2; x(2) + h(2)/2];
            x_ne = [x(1) + h(1)/2; x(2) + h(2)/2];
            % draw 4 edges
            plot([x_sw(1) x_se(1)],[x_sw(2) x_se(2)],'b-')
            hold on
            plot([x_nw(1) x_ne(1)],[x_nw(2) x_ne(2)],'b-')
            hold on
            plot([x_sw(1) x_nw(1)],[x_sw(2) x_nw(2)],'b-')
            hold on
            plot([x_se(1) x_ne(1)],[x_se(2) x_ne(2)],'b-')
            hold on
            % draw 1 facet
            fill([x_nw(1),x_sw(1),x_se(1),x_ne(1)]',[x_nw(2),x_sw(2),x_se(2),x_ne(2)]',color)
            hold on
        elseif length(N) == 3 || length(N) == 4
            % 3d cube plot
            x1 = [x(1)+h(1)/2; x(2)-h(2)/2; x(3)-h(3)/2];
            x2 = [x(1)+h(1)/2; x(2)+h(2)/2; x(3)-h(3)/2];
            x3 = [x(1)+h(1)/2; x(2)+h(2)/2; x(3)+h(3)/2];
            x4 = [x(1)+h(1)/2; x(2)-h(2)/2; x(3)+h(3)/2];
            x5 = [x(1)-h(1)/2; x(2)-h(2)/2; x(3)-h(3)/2];
            x6 = [x(1)-h(1)/2; x(2)+h(2)/2; x(3)-h(3)/2];
            x7 = [x(1)-h(1)/2; x(2)+h(2)/2; x(3)+h(3)/2];
            x8 = [x(1)-h(1)/2; x(2)-h(2)/2; x(3)+h(3)/2];
            % draw 6 facets
            fill3([x1(1),x2(1),x3(1),x4(1)]',[x1(2),x2(2),x3(2),x4(2)]',[x1(3),x2(3),x3(3),x4(3)]',color)
            hold on
            fill3([x5(1),x6(1),x7(1),x8(1)]',[x5(2),x6(2),x7(2),x8(2)]',[x5(3),x6(3),x7(3),x8(3)]',color)
            hold on
            fill3([x3(1),x2(1),x6(1),x7(1)]',[x3(2),x2(2),x6(2),x7(2)]',[x3(3),x2(3),x6(3),x7(3)]',color)
            hold on
            fill3([x4(1),x1(1),x5(1),x8(1)]',[x4(2),x1(2),x5(2),x8(2)]',[x4(3),x1(3),x5(3),x8(3)]',color)
            hold on
            fill3([x3(1),x4(1),x8(1),x7(1)]',[x3(2),x4(2),x8(2),x7(2)]',[x3(3),x4(3),x8(3),x7(3)]',color)
            hold on
            fill3([x2(1),x1(1),x5(1),x6(1)]',[x2(2),x1(2),x5(2),x6(2)]',[x2(3),x1(3),x5(3),x6(3)]',color)
        end
    end
    xc(index,:) = x';
end
%
if ~strcmp(edge,'EdgeOn')
    if length(N) == 2
        % project the plot onto x1-x2 space
        plot(xc(:,1),xc(:,2),'.',...
            'MarkerFaceColor',color,...
            'MarkerEdgeColor',color) % central points
    elseif length(N) == 3
        scatter3(xc(:,1),xc(:,2),xc(:,3),8,color,'filled','square')
    elseif length(N) == 4
        % 3d scatter to see all cells
        n = setdiff(1:4,[i,j,k]);
        scatter3(xc(:,i),xc(:,j),xc(:,k),18,xc(:,n),'filled','square')
    end
end
%
if length(N) == 2
    xlabel('$x_1$')
    ylabel('$x_2$')
    axis([lb(1) ub(1) lb(2) ub(2)]);
    box on
else
%     xlabel('$x_1$')
%     ylabel('$x_2$')
%     zlabel('$x_3$')
    
%     xlabel('$x(t)$')
%     ylabel('$y(t)$')
%     zlabel('$y(t-\tau)$')
    
    xlabel('$x(t)$')
    ylabel('$y(t)$')
    zlabel('$z(t)$')
    
    grid off
    axis([lb(i) ub(i) lb(j) ub(j) lb(k) ub(k)]);
    box on
end