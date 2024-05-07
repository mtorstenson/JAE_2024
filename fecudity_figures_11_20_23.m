%Creates figures used in manuscript

clear all; close all; clc

%Loads in data from a run of m_and_lh_sweep
load('3output_N=800_mx=0.8.mat','fecundity', 'migrant_m','one_stability_resident_popsize_mat', 'one_stability_migrant_popsize_mat','two_stability_resident_popsize_mat', 'two_stability_migrant_popsize_mat','delta_virulence_r_mat','proportion_migrant_susceptibles_mat','difference_in_eq_virulence_mat', 'one_stability_pop_virulence')
bfx = 1;
migrant_m = round(migrant_m,2);
delta_v_800_08 = delta_virulence_r_mat;
prop_m_800_08 = proportion_migrant_susceptibles_mat;
diff_800_08 = difference_in_eq_virulence_mat;
edv_800_08 = prop_m_800_08.*diff_800_08;
edv_n_800_08 = prop_m_800_08.*diff_800_08;
res_pop_diff_800_08 = two_stability_resident_popsize_mat - one_stability_resident_popsize_mat;
mig_pop_diff_800_08 = two_stability_migrant_popsize_mat - one_stability_migrant_popsize_mat;
orig_res_pop_800_08 = one_stability_resident_popsize_mat(1,1);
final_res_pop_800_08 = two_stability_resident_popsize_mat;
final_mig_pop_800_08 = two_stability_migrant_popsize_mat;
one_virulence_800_08 = difference_in_eq_virulence_mat + one_stability_pop_virulence;
one_stability_pop_virulence_800_08 = one_stability_pop_virulence;

tfsize = 12; % title size
fs1 = 10;  % axes labels
fs3 = 8; % axis numbering

%Creates figure for one population simulations
figure(1); clf
hh = gcf;
set(hh,'PaperUnits','centimeters');
set(hh,'Units','centimeters');
width = 7; height = 7;
xpos = 5;
ypos = 4;
set(gcf,'Position',[xpos ypos width height])

% size for individual panels
w = 0.75; % width
h = 0.75; % height
sy = 0.2; % starting x position
sx = 0.2; % starting y position

axes('position',[sx sy w h])
plot(fecundity,one_virulence_800_08(migrant_m == 0.8,:)-1,'-','Linewidth', 2)
yline(one_virulence_800_08(find(migrant_m==0.8),fecundity==1)-1,'-',{'m=0.8, f=1'})
ylim([0,10])
ylabel('Pathogen Strategy','FontSize',fs1)
xlabel('Fecundity Rate (f)','FontSize',fs1)
title('','FontSize',tfsize)
set(gca,'FontSize',fs3);

axes('position',[0 0 1 1],'visible','off')
    hold on;
    axis([0 1 0 1])
    
%creates figure where migrants and residents differ only in fecundity rate
figure(2); clf
hh = gcf;
set(hh,'PaperUnits','centimeters');
set(hh,'Units','centimeters');
width = 7; height = 12;
xpos = 5;
ypos = 4;
set(gcf,'Position',[xpos ypos width height])

% size for individual panels
w = 0.75; % width
h1 = 0.2; % height
h2 =0.5;222
d = 0.12; % spacing between
sy = 0.08; % starting x position
sx = 0.20; % starting y position

axes('position',[sx sy w h2])
hold on
plot(fecundity,final_res_pop_800_08(migrant_m==0.8,:)','-','Linewidth', 2)
plot(fecundity,final_mig_pop_800_08(migrant_m==0.8,:)','-','Linewidth', 2)
xline(bfx,'-',{'f_R'})
yline(0)
yline(orig_res_pop_800_08,'-','N_R')
ylabel('Population Size','FontSize',fs1)
xlabel('Fecundity rate of migrants','FontSize',fs1)
title('','FontSize',tfsize)
lgd = legend('Resident','Migrant','Location', 'southoutside');
title(lgd, 'Population')
set(gca,'FontSize',fs3);
box on
hold off

axes('position',[sx sy+h2+d w h1])
plot(fecundity,delta_v_800_08(migrant_m==0.8,:)','-','Linewidth', 2)
xline(bfx,'-',{'f_R'})
yline(0)
ylabel({'Change in','Pathogen Strategy'},'FontSize',fs1)
xlabel('Fecundity rate of migrants','FontSize',fs1)
title('','FontSize',tfsize)
set(gca,'FontSize',fs3);

% label subpanels
    axes('position',[0 0 1 1],'visible','off')
    hold on;
    text(0.02,0.14+h1+h2+d,'a)','FontSize',fs1);
    text(0.02,0.14+h2,'b)','FontSize',fs1);
axis([0 1 0 1])

% Backup previous settings
prePaperType = get(hh,'PaperType');
prePaperPosition = get(hh,'PaperPosition');
prePaperSize = get(hh,'PaperSize');

% Make changing paper type possible
set(hh,'PaperType','<custom>');
% Set the page size and position to match the figure's dimensions
position = get(hh,'Position');
set(hh,'PaperPosition',[0,0,position(3:4)]);
set(hh,'PaperSize',position(3:4));

print -dtiff -r600 3Fig1.tiff % save as 600dpi tiff
saveas(1,strcat(['3Fig1.jpg']))
print -dpdf -r600 3Fig1.pdf

print -dtiff -r600 3Fig2.tiff % save as 600dpi tiff
saveas(2,strcat(['3Fig2.jpg']))
print -dpdf -r600 3Fig2.pdf

% Restore the previous settings
set(hh,'PaperType',prePaperType);
set(hh,'PaperPosition',prePaperPosition);
set(hh,'PaperSize',prePaperSize);