function [N_base,d_base,mx_base,bfx_base,betamax_base,bmx_base,g_base,eps1_base,eps2_base] = set_baseline_parameter_values
% set default values to be used for all varied parameters

% set baseline parameter values
N_base = 800; %equilibrium population size w/o infection
d_base = 0.001; %density dependence parameter
mx_base = 0.8; %resident tolerance parameter
bfx_base = 1; %baseline fecundity rate of non-migratory population
betamax_base = 0.005; %maximum transmisson rate
bmx_base = 0.1; %resident baseline mortality
g_base = 0.1; % recovery rate for residents and migrants
eps1_base = 0.01; % use this for both em and ef
eps2_base = 0.001; % difference between migrant and resident mortality
