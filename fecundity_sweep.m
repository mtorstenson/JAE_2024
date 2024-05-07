%Sweeps through different migrant fecundity rates and outputs change in virulence
%for a 2 population model. The migratory population is added in after resident
%population reaches stability. 

N = 800; %equilibrium population size w/o infection
d = 0.001; %density dependence parameter
mx = 0.8; %resident tolerance parameter (index =17)
migrant_m = [0.8]; % migrant tolerance parameters 
fecundity = [0.1:0.1:2.0]; %migrant fecundity rates  
bmy = 0.09;
h = length(fecundity); 
count0 = 1;

%parameters for resident population (fixed throughout simulation)
betamax = 0.005; %maximum transmisson rate
beta = [0:(betamax/10):betamax]'; %transmission rates
kx = 0; %gamma dependence on beta
gammax = [.1+(kx*beta)]; %relatioship between transmission rate and probability of recovery (because kx = 0, gamma is the same for all values of beta)
lx = 0; %fecundity dependence on beta (fecundity the same for all values of beta because lx=0)
bfx = 0.8; %baseline fecundity rate of non-migratory population 
fecunditiesx = [bfx;[(bfx-(lx*beta))]]; %fecundity of susceptible and recovered classes
bmx = 0.09; %baseline mortality rate
miux = [bmx;[(sqrt(bmx)+(mx*(beta/betamax))).^(2)]]; %mortality rate of susceptible and infected classes, resident population, quadratic relationship between beta and miu

%parameters for migratory population
ky = kx; %gamma dependence on beta (same as in resident population)
gammay = [.1+(ky*beta)]; %relatioship between transmission rate and probability of recovery
ly = lx; %fecundity dependence on beta (same as in resident population)

%simulation parameters
p=0.001; %probability of mutation
years = 10000; %run will get cut off after 10,000 years if stability has not yet been reached
epsilon = 0.001; %stability cutoff

%initialize matrices of changes in virulence
delta_virulence_r_mat = nan(length(migrant_m),length(fecundity)); %for residents
delta_virulence_m_mat = nan(length(migrant_m),length(fecundity)); %for migrants

%initialize matrix of ratio of migrant:resident susceptibles at equilibrium
proportion_migrant_susceptibles_mat = nan(length(migrant_m),length(fecundity));

%initialize matrix of difference between strategy that evolves alone in migrants
%and residents
difference_in_eq_virulence_mat = nan(length(migrant_m),length(fecundity));

%initialize matrix of the weighted average of the optimal strategy
weighted_average_eq_virulence_mat = nan(length(migrant_m),length(fecundity));

%initialize matrix of resident population size at equilibrium when run alone with
%infection
one_stability_resident_popsize_mat = nan(length(migrant_m),length(fecundity));

%initialize matrix of migrant population size at equilibrium when run alone without
%infection
one_stability_migrant_popsize_mat = nan(length(migrant_m),length(fecundity));

%initialize matrix of resident population size at equilibrium when run with
%migrant and infection
two_stability_resident_popsize_mat = nan(length(migrant_m),length(fecundity));

%initialize matrix of migrant population size at equilibrium when run with
%resident and infection
two_stability_migrant_popsize_mat = nan(length(migrant_m),length(fecundity));

m_count = 1; %start count of migrant m values used

%run resident population to stabiltiy w/ pathogen

pop = [100;0;0;0;0;0;10;0;0;0;0;0]; % Initial population size for resident population
stability = 12; %resets stability parameter

%Run 1 population simulation of resident population with mutation
[stabilityPop, stabilityYear] = onePopSim_m(pop,beta,gammax,fecunditiesx,miux,p,years,epsilon,stability);
%tracks population and year when stability is reached
one_stability_pop = stabilityPop';
one_stability_year = stabilityYear;
[M,I] = max(stabilityPop);
one_stability_pop_max_value = M;
one_stability_pop_max_index = I;

pop1 = one_stability_pop'; % Initial population size for non-migratory population in this step is taken to be the end population from the one population simulation

%indexes through different values of migrant tolerance to infection
for my = migrant_m
    migratory_baseline_pop = nan(h,12); %initialize matrix to track migrant population before spillover (without pathogen)
    mig_stability_year = nan(h,1); %initialize matrix to track year at which stability reached (without pathogen)

    two_stability_pop_1 = nan(h,12); %initialize matrix to track equilibrium resident population after spillover
    two_stability_pop_1_max_value = nan(h,1); %initialize matrix to identify what infection class dominates in resident population after spillover
    two_stability_pop_1_max_index = nan(h,1); %initialize matrix to track what infection class dominates in resident population after spillover
    two_stability_pop_2 = nan(h,12); %initialize matrix to track equilibrium migrant population after spillover
    two_stability_pop_2_max_value = nan(h,1); %initialize matrix to identify what infection class dominates in migrant population after spillover
    two_stability_pop_2_max_index = nan(h,1); %initialize matrix to track what infection class dominates in migrant population after spillover
    two_stability_year = nan(h,1); %%initialize matrix to track when stability is reached in two-population simulation
    migrant_stability_pop = nan(h,12); %initialize matrix to track equilibrium migrant population when run alone with pathogen
    migrant_stability_year = nan(h,1); %initialize matrix year at which stability is reached in migrant population when run alone with pathogen

    lh_count = 1;     %Starts count to index life history speed values
    
    %indexes through different values of migrant life history speed
    for bfy = fecundity
        fecunditiesy = [bfy;[(bfy-(ly*beta))]]; %fecundity of migrant susceptible and infected classes

        %Run 'migratory' population to stability without pathogen
        miuy = [bmy;(sqrt(bmy)+(my*(beta/betamax))).^(2)]; %mortality rate for different classes of migrant
        pop = [110;0;0;0;0;0;0;0;0;0;0;0]; % Initial population size for migratory population (no infected individuals)
        stability = 12;
        [stabilityPop, stabilityYear] = onePopSim_m(pop,beta,gammay,fecunditiesy,miuy,p,years,epsilon,stability); %simulate migrant population without infection
        migratory_baseline_population = stabilityPop;

        migratory_baseline_pop(lh_count,:) = migratory_baseline_population';
        mig_stability_year(lh_count) = stabilityYear;

        %Add in migratory population, run without mutation

        %new simulation parameters
        pop2 = migratory_baseline_population; % Initial population of migratory population for two-population simulation (taken from equilibrium when run alone without infection)
        pop = vertcat(pop1, pop2); %combined resident and migrant population

        %re-set stability
        stability = 24;

        %Run two population simulation without mutation
        [stabilityPop, stabilityYear] = twoPopSim_m(pop,beta,gammax,gammay,fecunditiesx,fecunditiesy,miux,miuy,p,years,epsilon,stability);

        %tracks population and year when stability is reached
        two_stability_pop_1(lh_count,:) = stabilityPop(1:12)'; %resident
        two_stability_pop_2(lh_count,:) = stabilityPop(13:24)'; %migrant
        two_stability_year(lh_count) = stabilityYear; %year
        [M,I] = max(stabilityPop(1:12));
        two_stability_pop_1_max_value(lh_count) = M; %number of individuals in largest resident infection class
        two_stability_pop_1_max_index(lh_count) = I; %index largest resident infection class
        [M,I] = max(stabilityPop(13:24));
        two_stability_pop_2_max_value(lh_count) = M; %number of individuals in largest migrant infection class
        two_stability_pop_2_max_index(lh_count) = I; %index of largest migrant infection class
        
        %Checks what equilibrium virulence would be for migrant population alone
        pop = [100;0;0;0;0;0;10;0;0;0;0;0]; % Initial population size for migratory population
        stability = 12; %resets stability parameter
        [stabilityPop, stabilityYear] = onePopSim_m(pop,beta,gammay,fecunditiesy,miuy,p,years,epsilon,stability);
        %tracks population and year when stability is reached
        migrant_stability_pop(lh_count,:) = stabilityPop';
        migrant_stability_year(lh_count) = stabilityYear;

        lh_count = lh_count+1;

    end

    %Tracks population sizes at stability
    one_stability_popsize = sum(one_stability_pop,2)*ones(1,length(fecundity)); %resident alone
    two_stability_popsize_1 = sum(two_stability_pop_1,2); %resident together
    two_stability_popsize_2 = sum(two_stability_pop_2,2); %migrant together

    %Determines average pathogen strategy
    one_stability_pop_virulence = (one_stability_pop(2:end)*(1:11)')./sum(one_stability_pop(2:end),2); %resident alone
    migrant_stability_virulence = (migrant_stability_pop(:,2:end)*(1:11)')./sum(migrant_stability_pop(:,2:end),2); %migrant alone
    two_stability_pop_1_virulence = (two_stability_pop_1(:,2:end)*(1:11)')./sum(two_stability_pop_1(:,2:end),2); %resident with migrant
    two_stability_pop_2_virulence = (two_stability_pop_2(:,2:end)*(1:11)')./sum(two_stability_pop_2(:,2:end),2); %migrant with resident
    delta_virulence_r = two_stability_pop_1_virulence - one_stability_pop_virulence; %change in resident virulence
    delta_virulence_m = two_stability_pop_2_virulence - one_stability_pop_virulence; %difference between migrant virulence and resident alone virulence
    difference_virulence = migrant_stability_virulence - one_stability_pop_virulence; %difference between migrant alone and resident alone virulence
    
    migrant_proportion = two_stability_pop_2(:,1)./(two_stability_pop_1(:,1)+two_stability_pop_2(:,1)); %identifies proportion of susceptibles individuals that are migrants
    
    %puts values from this set of parameters into matrices for the entire sweep
    delta_virulence_r_mat(m_count,:) = delta_virulence_r';
    delta_virulence_m_mat(m_count,:) = delta_virulence_m';
    difference_in_eq_virulence_mat(m_count,:) = difference_virulence';
    proportion_migrant_susceptibles_mat(m_count,:) = migrant_proportion;
    weighted_average_eq_virulence_mat(m_count,:) = migrant_proportion' .* migrant_stability_virulence' + (1-migrant_proportion)' * one_stability_pop_virulence';
    one_stability_resident_popsize_mat(m_count,:) = one_stability_popsize;
    one_stability_migrant_popsize_mat(m_count,:) = sum(migratory_baseline_pop,2);
    two_stability_resident_popsize_mat(m_count,:) = two_stability_popsize_1;
    two_stability_migrant_popsize_mat(m_count,:) = two_stability_popsize_2;

    m_count = m_count+1;
end

%Saves the output of the script as a .mat file
save output.mat
save(strcat(['3output_N=' num2str(N) '_mx=' num2str(mx)  '.mat']))