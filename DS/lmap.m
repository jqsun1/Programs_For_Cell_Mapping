%Logistic Map
swEPSfigure;


L = @(x,u) u*x*(1-x);
close all
u = 2;
x = 0:0.01:1;
for i=1:length(x)
   Lx(i) = L(x(i), u); 
end
figure(1)
plot(x,Lx,'k');
hold on
plot(x,x,'--k');

set(gca,'XTick',[0 0.5])
set(gca,'XTickLabel', {'$x_1^*$', '$x_2^*$'}, 'TickLabelInterpreter','latex')
set(gca,'defaulttextinterpreter','latex');
set(gcf,'defaulttextinterpreter','latex');
xlabel(strcat('$x$'),'fontsize',18);
ylabel(strcat('$f$ and $g$'),'fontsize',18);
% set(gcf, 'Units','centimeters', 'Position',[1 1 9 9]);
% set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[18, 18]);
set(gcf,'color','white')
print -depsc u2.eps

u = 3.4;
x = 0:0.01:1;
for i=1:length(x)
   Lx(i) = L(x(i), u); 
end
figure(2)
plot(x,Lx,'k');
hold on
plot(x,x,'--k');

x1 = 0;
x2 = 0.7058823529411737;
x3 = 0.4519632476261533;
x4 = 0.842154399432673;
fx3= L(x3,u);
fx4= L(x4,u);

v = [x3 x4  x4 x3 x3];
fv= [x3 fx4 x4 fx3 x3];
plot(v,fv,'-.k')

set(gca,'XTick',[x1 x3 x2 x4])
set(gca,'XTickLabel', {'$x_1^*$', '$x_2^*$','$x_3^*$', '$x_4^*$'}, 'TickLabelInterpreter','latex')

% plot(x3,fx3,'*')
% plot(x4,fx4,'*')
% plot(x3,x3,'*')
% plot(x4,x4,'*')

xlabel(strcat('$x$'),'fontsize',18);
ylabel(strcat('$f$ and $g$'),'fontsize',18);
% set(gcf, 'Units','centimeters', 'Position',[1 1 9 9]);
% set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[18, 18]);    
set(gcf,'color','white')
print -depsc u3.eps

