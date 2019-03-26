function plot_beam_shape(xf_lb,xf_ub,Lb,v)
% -------------------------------------------------------------------------
% Plot the step beams with different geometric configurations.
%
% Input arguments:
%       Lb:    length of beam
%       He:    each row represents the different heights of each segment
%
% By: Free Xiong; 2014-11-25
% -------------------------------------------------------------------------
swEPSfigure;
figure
Nseg = 50;
N_spline = length(v);
xs = ((1:N_spline)'-1)*Lb/(N_spline-1);
ys = v;
xt = Lb/Nseg/2 + ((1:Nseg)'-1)*Lb/Nseg;
He(1,:) = spline(xs,v(1,:),xt);
He(2,:) = spline(xs,v(2,:),xt);
[m,n] = size(He);
for i = 1:m
    subplot(m,1,i)
    for j = 1:n
        hj = He(i,j);
        plot([(j-1)*Lb/n, j*Lb/n],[-hj/2, -hj/2],'b');
        hold on
        plot([(j-1)*Lb/n, j*Lb/n],[hj/2, hj/2],'b');
        hold on
        plot([(j-1)*Lb/n, (j-1)*Lb/n],[-hj/2, hj/2],'b');
        hold on
        plot([j*Lb/n, j*Lb/n],[-hj/2, hj/2],'b');
        hold on
        ylim = get(gca,'Ylim');
        set(gca,'YLim',[-0.02, 0.02]);
        plot([xf_lb xf_lb],get(gca,'Ylim'),'r')
        plot([xf_ub xf_ub],get(gca,'Ylim'),'r')
%         set(gca,'xtick',[])
%         set(gca,'xticklabel',[])
        set(gca,'ytick',[])
        set(gca,'yticklabel',[])
        title(['Beam No. ',num2str(i)])
    end
end
box on