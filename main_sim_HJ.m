% ========================================================================
% MONTE CARLO SIMULATION
% SECTION 3.1: AR(1) MODEL
% ------------------------------------------------------------------------
% Appendix B.1 (include VAR(1) in pool):
%   spec_fn = @() spec_HJ_VAR;
% Appendix B.2 (omitted variable bias):
%   sim_par.s_eu = -0.05 or 0.05;
%   sim_fn = @() sim_HJ_OV(sim_par);
% ========================================================================

clear
clc
addpath(genpath(pwd))

% FILE NAME TO SAVE OUTPUT
out_name = 'sim_HJ';

% SIMULATION PARAMETERS
sim_par.rho = 0.95;     % persistance of AR(1) and std dev of shocks
sim_par.T = 150;    	% number of periods
n_sim = 5e5;            % number of simulations

% FUNCTIONS
sim_fn  = @() sim_HJ(sim_par);          % simulation
spec_fn = @() spec_HJ;                  % model and estimation specification
sel_fn  = @(est, wt, spec)...
    sel_out_sim_HJ(est, wt, spec);      % select output to save

% VARIABLES OF INTEREST
var_list = {'y'};

% MAIN COMPUTATIONS
calc_simulation;