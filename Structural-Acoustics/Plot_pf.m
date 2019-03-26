clear
clc
close all
fn = 'new_fem_pop120_gen120_sound_omega_SS_250_450_';
% fn = 'beam_data\fem_pop100_gen80_vs_CC_TMM_';
ps_tot = [];
pf_tot = [];
time_avg = 0;
for i = 1:10
    filename = [fn,num2str(i)];
    load(filename);
    ps_tot = [ps_tot; ps_PSO];
    pf_tot = [pf_tot; pf_PSO];
    time_avg = time_avg + time_MOPSO;
end
time_avg = time_avg/10;

% ps_tot = ps_tot(pf_tot(:,2)<2e-4,:);
% pf_tot = pf_tot(pf_tot(:,2)<2e-4,:);
[pf_all, ps_all] = dom_chk(pf_tot,ps_tot);
% load('new_fem_pop120_gen120_sound_omega_SS_250_450_10.mat')



figure
scatter(10*log(pf_all(:,1)/200/1e-12)/log(10),pf_all(:,2)/2/pi,12,'filled');
hold on
plot([0 100],[-19.62 -19.62],'r');
hold on
% plot([1.4e-3 1.4e-3],[-75 -40],'r');
% scatter(pf_all(:,1),-pf_all(:,3)/2/pi,12,'filled');
axis([46 64 -22 -18])
set(gca,'xtick',[46:2:64])
box on
xlabel('Radiated Sound Power (dB)')
ylabel('$-\omega_1$ (Hz)')
swFigSize
swEPSfigure
print('-depsc','figures\SS_250_450_Fundamental_Frequency_opt_hand_pareto_front')

%%
clear
fn = 'new_fem_pop120_gen120_sound_omega_CC_250_450_';
% fn = 'beam_data\fem_pop100_gen80_vs_CC_TMM_';
ps_tot = [];
pf_tot = [];
time_avg = 0;
for i = 1:10
    filename = [fn,num2str(i)];
    load(filename);
    ps_tot = [ps_tot; ps_PSO];
    pf_tot = [pf_tot; pf_PSO];
    time_avg = time_avg + time_MOPSO;
end
time_avg = time_avg/10;

% ps_tot = ps_tot(pf_tot(:,2)<2e-4,:);
% pf_tot = pf_tot(pf_tot(:,2)<2e-4,:);
[pf_all, ps_all] = dom_chk(pf_tot,ps_tot);
% load('new_fem_pop120_gen120_sound_omega_SS_250_450_10.mat')


figure
scatter(10*log(pf_all(:,1)/200/1e-12)/log(10),pf_all(:,2)/2/pi,12,'filled');
hold on
plot([0 100],[-44.46 -44.46],'r');
hold on
% plot([1.4e-3 1.4e-3],[-75 -40],'r');
% scatter(pf_all(:,1),-pf_all(:,3)/2/pi,12,'filled');
 axis([46 58 -75 -40])
set(gca,'xtick',[46:2:64])
box on
xlabel('Radiated Sound Power (dB)')
ylabel('$-\omega_1$ (Hz)')

swFigSize
swEPSfigure
print('-depsc','figures\CC_250_450_Fundamental_Frequency_opt_hand_pareto_front')

%%
clear
load('Pso_scm_nonunifrom_beam_input_10_div10x10_fem_maxiter_10.mat')
% close all
% scatter(pf_all(:,1),pf_all(:,2),12,'filled');
% hold on
% scatter(pf(:,1),pf(:,2),12,'filled');

figure
scatter(pf(:,1),10*log(pf(:,2)/1e-12)/log(10),30,'o')
hold on
scatter(pf_all(:,1),10*log(400*pf_all(:,2)/1e-12)/log(10),30,'filled')
ylabel('Radiated Sound Power (dB)')
xlabel('Mass (kg)')
box on
swFigSize
swEPSfigure
print('-depsc','figures\pso_scm_with_MOPSO')
%%


% 
clear
% close all
figure


load('fem_pop100_gen80_vs_CC_TMM_1.mat');
plot(pf_PSO(:,1),10*log(pf_PSO(:,2)/1e-12)/log(10),'bo','MarkerFaceColor','b','markersize',6);
hold on
load('CC_pop100_gen80_freq_100_to_300_jizhong_all.mat')
plot(pf_all(:,1),10*log(pf_all(:,2)/200/1e-12)/log(10),'ro','MarkerFaceColor','w','markersize',6);
xlabel('Mass (kg)')
ylabel('Radiated Sound Power (dB)')
axis([18 70 50 80])
swFigSize
swEPSfigure
print('-depsc','figures\FEM_vs_TMM_CC')

%%
clear
figure
swEPSfigure
ps_tot = [];
pf_tot = [];
for kkk = 1 :8
    
    filename='fem_pop100_gen80_vs_TMM_';
    s1 = num2str(kkk);
    filename1 = [filename,s1];
    load(filename1);    
    ps_tot = [ps_tot; ps_PSO];
    pf_tot = [pf_tot; pf_PSO];

end
[pf, ps] = dom_chk(pf_tot,ps_tot);
% load('beam_data\fem_pop100_gen80_vs_TMM_8.mat');
plot(pf(:,1),10*log(pf(:,2)/1e-12)/log(10),'bo','MarkerFaceColor','b','markersize',6);
hold on
load('SS_pop100_gen80_freq_200_to_600_jizhong_all.mat')
plot(pf_all(:,1),10*log(pf_all(:,2)/400/1e-12)/log(10),'ro','MarkerFaceColor','w','markersize',6);
xlabel('Mass (kg)')
ylabel('Radiated Sound Power (dB)')
axis([33 60 55 70])
swFigSize


print('-depsc','figures\FEM_vs_TMM_SS')

%%
clear
clc
close all
% load('Pso_scm_nonunifrom_beam_input_8_div10x10_fem_maxiter_10')
% v = [ps(11,:);ps_PSO(15,:)];
load('Pso_scm_nonunifrom_beam_input_10_div10x10_fem_maxiter_10')

v = [ps(1,:);ps(end,:)];
% 88 dB 89.58----91.13  92.97.....35.23 35.13

plot_beam_shape(0.7071*Lb,0.7071*Lb,Lb,v)

% saveas(gcf,'figures\case__psoscm_vs_pso.eps')
% print('-depsc','figures\case__psoscm_vs_pso_1')
% print('-depsc','figures\extreme_design_SS_for_sound_and_mass_pso_scm_f_200_600')
