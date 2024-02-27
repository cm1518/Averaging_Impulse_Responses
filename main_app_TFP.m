% ========================================================================
% EMPIRICAL APPLICATION
% SECTION 4.2: TFP SHOCK
% ========================================================================

clear
clc
addpath(genpath(pwd))

% READ DATA
data = data_TFP;

% MODEL AND ESTIMATION SPECIFICATIONS
spec_TFP;

% VARIABLES TO COMPUTE WEIGHTS FOR
var_list = {'unemp'};

% FILE NAME TO SAVE OUTPUT
out_name = 'TFP';

% MAIN COMPUTATIONS
calc_application;
