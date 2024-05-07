clear all; clc
% run sensitivity analysis using Latin Hypercube Sampling (LHS) and PRCC
% This code relies on previously published code:
%   it uses the script LHS_Call.m available from Simeone Marino here:
%   http://malthus.micro.med.umich.edu/lab/usadata/data/matlab-ode-lhs.10-03-2023.zip
%   more information is available from this webpage:
%   http://malthus.micro.med.umich.edu/lab/usanalysis.html
%   and this publication:
%   Marino S, Hogue IB, Ray CJ, Kirschner DE (2008) A methodology for performing 
%        global uncertainty and sensitivity analysis in systems biology. J Theor
%        Biol 254(1):178–196
%        https://www.sciencedirect.com/science/article/pii/S0022519308001896
%
% output measure: change in pathogen virulence infecting a resident
%   population in isolation vs in the presence of a migratory population
%   (delta_virulence_R)
%
% we consider four sets of sensitivity analysis runs:
%   1. set em=0 and vary ef < 0 and emiu < 0
%   2. set em=0 and vary ef > 0 and emiu > 0
%   3. set ef=0 and emiu = 0 and vary em < 0
%   4. set ef=0 and emiu = 0 and vary em > 0

runs = 10000; %Number of samples to run

% set default values to be used for all varied parameters
[N_base,d_base,mx_base,bfx_base,betamax_base,bmx_base,g_base,eps1_base,eps2_base] = set_baseline_parameter_values;

% use Latin Hypercube Sampling to generate parameter combinations
% (avoid setting minimum values as zero, otherwise distributions change)
%                 min-value base-value  max-value  #samples  distribution-type
v_N       = LHS_Call(400,    N_base,       1600, 0, runs, 'unif');
v_d       = LHS_Call(0.0001, d_base,       0.01, 0, runs, 'unif');
v_mx      = LHS_Call(0.01,   mx_base,      2,    0, runs, 'unif');
v_bfx     = LHS_Call(0.1,    bfx_base,     2,    0, runs, 'unif');
v_betamax = LHS_Call(0.003,  betamax_base, 0.01, 0, runs, 'unif');
v_bmx     = LHS_Call(0.01,   bmx_base,     0.2,  0, runs, 'unif');
v_g       = LHS_Call(0.05,   g_base,       0.2,  0, runs, 'unif');
v_eps1    = LHS_Call(0.001,  eps1_base,    0.1,  0, runs, 'unif');
v_eps2    = LHS_Call(0.0001, eps2_base,    0.01, 0, runs, 'unif');

% set up LHS matrix - this holds all the combinations of varied parameters
LHSmatrix = [v_N v_d v_mx v_bfx v_betamax v_bmx v_g v_eps1 v_eps2];

% set up for 4 sensitivity sets to run
deltavir1 = NaN(1,runs);
deltavir2 = NaN(1,runs);
deltavir3 = NaN(1,runs);
deltavir4 = NaN(1,runs);

tic
for x = 1:runs
    
    x
    
    % pull the appropriate set of parameters
    N       = LHSmatrix(x,1); %equilibrium population size w/o infection
    d       = LHSmatrix(x,2); %density dependence parameter
    mx      = LHSmatrix(x,3); %resident tolerance parameter
    bfx     = LHSmatrix(x,4); %baseline fecundity rate of non-migratory population
    betamax = LHSmatrix(x,5); %maximum transmisson rate
    bmx     = LHSmatrix(x,6); %resident baseline mortality
    g       = LHSmatrix(x,7); % recovery rate for residents and migrants
    
    
    %   1. set em=0 and vary ef < 0 and emiu < 0
    em   = 0; % difference between migrant and resident tolerance
    ef   = -LHSmatrix(x,8); % difference between migrant and resident fecundity
    emiu = -LHSmatrix(x,9); % difference between migrant and resident mortality
    [delta_virulence_r] = run_single_simulation(N,d,mx,em,bfx,ef,betamax,bmx,emiu,g); % run simulation
    deltavir1(x) = delta_virulence_r; % save output measure

    %   2. set em=0 and vary ef > 0 and emiu > 0
    em   = 0; % difference between migrant and resident tolerance
    ef   = LHSmatrix(x,8); % difference between migrant and resident fecundity
    emiu = LHSmatrix(x,9); % difference between migrant and resident mortality
    [delta_virulence_r] = run_single_simulation(N,d,mx,em,bfx,ef,betamax,bmx,emiu,g); % run simulation
    deltavir2(x) = delta_virulence_r; % save output measure

    %   3. set ef=0 and emiu = 0 and vary em < 0
    em   = -LHSmatrix(x,8); % difference between migrant and resident tolerance
    ef   = 0; % difference between migrant and resident fecundity
    emiu = 0; % difference between migrant and resident mortality
    [delta_virulence_r] = run_single_simulation(N,d,mx,em,bfx,ef,betamax,bmx,emiu,g); % run simulation
    deltavir3(x) = delta_virulence_r; % save output measure
    
    %   4. set ef=0 and emiu = 0 and vary em > 0
    em   = LHSmatrix(x,8); % difference between migrant and resident tolerance
    ef   = 0; % difference between migrant and resident fecundity
    emiu = 0; % difference between migrant and resident mortality
    [delta_virulence_r] = run_single_simulation(N,d,mx,em,bfx,ef,betamax,bmx,emiu,g); % run simulation
    deltavir4(x) = delta_virulence_r; % save output measure

save sensitivity.mat *
    
end
toc

%%
load sensitivity.mat

per1 = round(sum(deltavir1<0)/runs,2)*100;
per2 = round(sum(deltavir2>0)/runs,2)*100; 
per3 = round(sum(deltavir3>0)/runs,2)*100;
per4 = round(sum(deltavir4<0)/runs,2)*100;

% set plotting parameters
fs1 = 10;  % axes labels
fs2 = 9;  % axis numbering
lw2 = 1; % fig edges

width = 14;
height = 14;
xpos = 3;
ypos = 2;
sx = 0.12;
sy = 0.08;
w = 0.35;
he = 0.32;
dx = 0.16;
dy = 0.18;


figure(1); clf
hh = gcf;
set(hh,'PaperUnits','centimeters');
set(hh,'Units','centimeters');
set(gcf,'Position',[xpos ypos width height])



axes('position',[sx sy+dy+he w he])
    histogram(deltavir1);
    title({'em=0, ef < 0, emiu < 0',strcat(['(' num2str(per1) '% negative)'])})
    xlabel('Change in virulence','FontSize',fs1);
    ylabel('Number of simulations','FontSize',fs1)

axes('position',[sx+dx+w sy+dy+he w he])
    histogram(deltavir2);
    title({'em=0, ef > 0, emiu > 0',strcat(['(' num2str(per2) '% positive)'])})
    xlabel('Change in virulence','FontSize',fs1);
    ylabel('Number of simulations','FontSize',fs1)

axes('position',[sx sy w he])
    histogram(deltavir3);
    title({'ef=0, emiu=0, em < 0 ',strcat(['(' num2str(per3) '% positive)'])})
    xlabel('Change in virulence','FontSize',fs1);
    ylabel('Number of simulations','FontSize',fs1)

axes('position',[sx+dx+w sy w he])
    histogram(deltavir4); 
    title({'ef=0, emiu=0, em > 0 ',strcat(['(' num2str(per4) '% negative)'])})
    xlabel('Change in virulence','FontSize',fs1);
    ylabel('Number of simulations','FontSize',fs1)
    set(gca,'FontSize',fs2,'LineWidth',lw2,'Fontname', 'Arial'); 

% label subpanels
axes('position',[0 0 1 1],'visible','off')
hold on
     text(0.02,0.98,'a)','horizontalalignment','center')
     text(0.52,0.98,'b)','horizontalalignment','center')
     text(0.02,0.46,'c)','horizontalalignment','center')
     text(0.52,0.46,'d)','horizontalalignment','center')
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

print -djpeg -r600 FigureS6_sensitivity.jpg

% Restore the previous settings
set(hh,'PaperType',prePaperType);
set(hh,'PaperPosition',prePaperPosition);
set(hh,'PaperSize',prePaperSize);



