% ========================================================================
% MONTE CARLO SIMULATION
% APPENDIX C: NEW KEYNESIAN MODEL FROM SMETS AND WOUTERS (2007)
% ========================================================================

clear
clc
addpath(genpath(pwd))

% FILE NAME TO SAVE OUTPUT
out_name = 'sim_SW';

% SIMULATION PARAMETERS
sim_par.T_sim = 150;    % number of periods for main sample
sim_par.T_pre = 100;    % number of periods for presample
n_sim = 2.5e4;          % number of simulations

% FUNCTIONS
sim_fn  = @() sim_SW(sim_par);                              % simulation
spec_fn = @() spec_SW;                                      % model and estimation specification
sel_fn  = @(est, wt, spec)sel_out_sim_SW(est, wt, spec);    % select output to save

% VARIABLES OF INTEREST
var_list = {'GDP'  'Inflation'  'FFR'  'Consumption'  'Investment'  'Wage'  'Hours'};

% MAIN COMPUTATIONS
calc_simulation;

