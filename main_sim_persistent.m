% ========================================================================
% MONTE CARLO SIMULATION
% SECTION 3.2: MODEL WITH AR(1) SHOCK
% ========================================================================

clear
clc
addpath(genpath(pwd))

% FILE NAME TO SAVE OUTPUT
out_name = 'sim_peristent';

% SIMULATION PARAMETERS
sim_par.rho_y = 0.97;   % persistence of y
sim_par.rho_u = 0.75;   % persistence of shock
sim_par.T_sim = 250;    % number of periods for main sample
sim_par.T_pre = 1e3;    % number of periods for presample
n_sim = 2.5e4;          % number of simulations

% FUNCTIONS
sim_fn  = @() sim_peristent(sim_par);               % simulation
spec_fn = @() spec_peristent;                       % model and estimation specification
sel_fn  = @(est, wt, spec)...
    sel_out_sim_peristent(est, wt, spec, sim_par);  % select output to save

% VARIABLES OF INTEREST
var_list = {'y'};

% MAIN COMPUTATIONS
calc_simulation;
