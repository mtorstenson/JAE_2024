function [delta_virulence_r] = run_single_simulation(N,d,mx,em,bfx,ef,betamax,bmx,emiu,g)
%Runs one simulation in which an equilibrium virulence is found for the
%resident population, then a migrant population is added in and the change
%in equilibrium virulence is identified 

gx = g; %resident recovery rate
gy = g; %migrant recovery rate


my = max(0,mx+em); % migrant tolerance parameters 
bfy = max(0,bfx+ef); %migrant fecundity rates 
bmy = max(0,bmx+emiu); 

%simulation parameters
p=0.001; %probability of mutation
years = 10000; %run will get cut off after 10,000 years if stability has not yet been reached
epsilon = 0.001; %stability cutoff

%set vectors
beta = [0:(betamax/10):betamax]'; %transmission rates
kx = 0; %gamma dependence on beta
ky = kx; %gamma dependence on beta (same as in resident population)
gammax = [gx+(kx*beta)]; %relatioship between transmission rate and probability of recovery (because kx = 0, gamma is the same for all values of beta)
gammay = [gy+(ky*beta)]; %relatioship between transmission rate and probability of recovery
lx = 0; %fecundity dependence on beta (fecundity the same for all values of beta because lx=0)
ly = 0; %fecundity dependence on beta (same as in resident population)
fecunditiesx = [bfx;[(bfx-(lx*beta))]]; %fecundity of susceptible and recovered classes
fecunditiesy = [bfy;[(bfy-(ly*beta))]]; %fecundity of migrant susceptible and infected classe
miux = [bmx;[(sqrt(bmx)+(mx*(beta/betamax))).^(2)]]; %mortality rate of susceptible and infected classes, resident population, quadratic relationship between beta and miu
miuy = [bmy;(sqrt(bmy)+(my*(beta/betamax))).^(2)]; %mortality rate for different classes of migrant

%run resident population to stabiltiy w/ pathogen
pop = [100;0;0;0;0;0;10;0;0;0;0;0]; % Initial population size for resident population
stability = 12; %resets stability parameter

%Run 1 population simulation of resident population with mutation
[stabilityPop, stabilityYear] = onePopSim_m(pop,beta,gammax,fecunditiesx,miux,p,years,epsilon,stability);
%tracks population and year when stability is reached
one_stability_pop = stabilityPop';

pop1 = one_stability_pop'; % Initial population size for non-migratory population in this step is taken to be the end population from the one population simulation

%Run 'migratory' population to stability without pathogen
pop = [110;0;0;0;0;0;0;0;0;0;0;0]; % Initial population size for migratory population (no infected individuals)
stability = 12;
[stabilityPop, stabilityYear] = onePopSim_m(pop,beta,gammay,fecunditiesy,miuy,p,years,epsilon,stability); %simulate migrant population without infection
pop2 = stabilityPop;

%Add in migratory population
pop = vertcat(pop1, pop2); %combined resident and migrant population

%re-set stability
stability = 24;

%Run two population simulation without mutation
[stabilityPop, stabilityYear] = twoPopSim_m(pop,beta,gammax,gammay,fecunditiesx,fecunditiesy,miux,miuy,p,years,epsilon,stability);

%tracks population and year when stability is reached
two_stability_pop_1 = stabilityPop(1:12)'; %resident
two_stability_pop_2 = stabilityPop(13:24)'; %migrant

%Determines average pathogen strategy
one_stability_pop_virulence = (one_stability_pop(2:end)*(1:11)')./sum(one_stability_pop(2:end),2); %resident alone
two_stability_pop_1_virulence = (two_stability_pop_1(:,2:end)*(1:11)')./sum(two_stability_pop_1(:,2:end),2); %resident with migrant
two_stability_pop_2_virulence = (two_stability_pop_2(:,2:end)*(1:11)')./sum(two_stability_pop_2(:,2:end),2); %migrant with resident
delta_virulence_r = two_stability_pop_1_virulence - one_stability_pop_virulence; %change in resident virulence
delta_virulence_m = two_stability_pop_2_virulence - one_stability_pop_virulence; %difference between migrant virulence and resident alone virulence
