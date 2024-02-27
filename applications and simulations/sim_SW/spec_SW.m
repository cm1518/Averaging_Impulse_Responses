function [spec, est] = spec_SW

% ========================================================================
% SELECT SPECIFICATIONS AND MODELS TO BE ESTIMATED
% >> APPENDIX C: MONTE CARLO WITH SMETS AND WOUTERS (2007) NK MODEL
% ========================================================================

% ------------------------------------------------------------------------
% MODELS, IDENTIFICATION, VARIABLES
% ------------------------------------------------------------------------

% MODELS
% options: 'VAR', 'LP', 'FAVAR', 'single equation'
spec.models = {'VAR', 'LP'};

% IDENTIFICATION SCHEMES
% >> VAR: 'Cholesky', 'internal IV'
% >> LP (frequentist): 'recursive', 'non-recursive'
% >> single equation: N/A, leave as ''
spec.id = {'internal IV', 'internal IV'};

% NUMBER OF LAGS, MODELS
spec.n_lag = [1,1];
spec.n_model = size(spec.models,2); 

% SELECT VARIABLES FOR EACH MODEL
est = cell(spec.n_model,1);
% >> 1. VAR
est{1}.var = {'GDP','Inflation','FFR','Consumption','Investment','Wage','Hours'};
est{1}.IV  = {'Instrument'};
% >> 2. LP
est{2}.var = est{1}.var;
est{2}.IV  = est{1}.IV;


% ------------------------------------------------------------------------
% ESTIMATION DETAILS
% ------------------------------------------------------------------------

% ESTIMATION PARAMETERS
spec.H        = 20;             % number of horizons
spec.n_draw   = 1e3;            % number of draws
spec.n_blk    = 2;              % number of blocks for sample splitting. 0 indicates you train and test on the entire dataset 
spec.constant = 'true';         % options: 'true', 'false'
spec.btstrp   = 'wild';         % options: 'standard', 'wild'; only matters for frequentist estimation
spec.uncon    = 'false';        % options: 'true', 'false'
spec.paradigm = 'Bayesian';     % options: 'Bayesian', 'frequentist'

for mm = 1:length(est)

    % LOAD MODEL SPECIFICATION
    est{mm}.model  = spec.models(mm);     	% model
    est{mm}.id     = spec.id(mm);           % identification
    est{mm}.n_lag  = spec.n_lag(mm);        % number of lags

    % VARIABLE TO NORMALIZE IRF TO
    % >> matters only for plotting IRF
    % >>  -1: 1 std dev; only for internal IV VAR;
    % >>   0: 1 unit change in instrument
    % >> i>0: 1 unit change in endogenous variable i
    est{mm}.norm_i = find(strcmp(est{mm}.var,'FFR')); % 1 unit change in FFR

end


end