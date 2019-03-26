clear all
close all

dtlz1k2 = [15.2755 31.8103 47.9015];
dtlz2k2 = [15.6557 28.8584 37.0169];
dtlz3k2 = [15.7142 30.4780 42.5710];
dtlz4k2 = [15.2500 31.4285 44.4541];
dtlz5k2 = [14.3442 28.4792 38.7896];
dtlz6k2 = [21.5322 60.0414 96.1908];
dtlz7k2 = [10.7633 11.1688 11.25];

dtlz1k3 = [14.8251 27.6752 40.7152];
dtlz2k3 = [14.3750 24.3462 33.1129];
dtlz3k3 = [15.4729 26.6037 35.9134];
dtlz4k3 = [17.0802 36.5410 52.3000];
dtlz5k3 = [14.2142 24.1197 32.5384];
dtlz6k3 = [18.9655 42.2727 70.1735];
dtlz7k3 = [12.2463 12.2169 10.2234];

x = [10 100 1000];

swEPSfigure

figure;
plot(x,dtlz1k2,'k-',x,dtlz2k2,'k--',x,dtlz3k2,'k-.',x,dtlz4k2,'k:',x,dtlz5k2,'b-',x,dtlz6k2,'b--',x,dtlz7k2,'b-.')
xlabel('$m$')
ylabel('Acceleration')
legend('DTLZ1','DTLZ2','DTLZ3','DTLZ4','DTLZ5','DTLZ6','DTLZ7')
set(gcf,'color','white')
pause
print -depsc acc-k2.eps


figure;
plot(x,dtlz1k3,'k-',x,dtlz2k3,'k--',x,dtlz3k3,'k-.',x,dtlz4k3,'k:',x,dtlz5k3,'b-',x,dtlz6k3,'b--',x,dtlz7k3,'b-.')
xlabel('$m$')
ylabel('Acceleration')
legend('DTLZ1','DTLZ2','DTLZ3','DTLZ4','DTLZ5','DTLZ6','DTLZ7')
set(gcf,'color','white')
pause
print -depsc acc-k3.eps

k = [2,3,4,5];
tp = [.1 .3 .5 .8]*10e5;
ts = [1.1 3 3.5 4.05]*10e5;

figure;
plot(k, tp,'k-', k, ts,'k--')
legend('Parallel','Sequential')
xlabel('Objectives')
ylabel('Time')
set(gcf,'color','white')
pause
print -depsc obj.eps

figure
ac = [12.8,12,8.2,5.6];
plot(k, ac,'k-')
xlabel('Objectives')
ylabel('Acceleration')
set(gcf,'color','white')
print -depsc accobj.eps

