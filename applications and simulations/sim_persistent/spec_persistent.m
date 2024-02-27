function [spec, est] = spec_persistent

% ========================================================================
% SELECT SPECIFICATIONS AND MODELS TO BE ESTIMATED
% >> SECTION 3.2: MONTE CARLO WITH PERSISTENT SHOCK
% ========================================================================

% ------------------------------------------------------------------------
% MODELS, IDENTIFICATION, VARIABLES
% ------------------------------------------------------------------------

% MODELS
% options: 'VAR', 'LP', 'FAVAR', 'single equation'
spec.models = {'VAR', 'LP'};

% IDENTIFICATION SCHEMES
% >> VAR: 'Cholesky', 'internal IV'
% >> LP: 'recursive', 'non-recursive'
% >> single equation: N/A, leave as ''
spec.id = {'internal IV', 'non-recursive'};

% NUMBER OF LAGS, MODELS
spec.n_lag = [1,1];
spec.n_model = size(spec.models,2); 

% SELECT VARIABLES FOR EACH MODEL
est = cell(spec.n_model,1);
% >> 1. VAR: internal instrument
est{1}.var = {'y'};
est{1}.IV  = {'v'};
% >> 2. LP: non-recursive
est{2}.var = {'y'};
est{2}.IV  = {'v'};


% ------------------------------------------------------------------------
% ESTIMATION DETAILS
% ------------------------------------------------------------------------

% ESTIMATION PARAMETERS
spec.H        = 20;             % number of horizons
spec.n_draw   = 1e3;            % number of draws
spec.n_blk    = 2;              % number of blocks for sample splitting. 0 indicates you train and test on the entire dataset 
spec.constant = 'false';        % options: 'true', 'false'
spec.btstrp   = 'wild';         % options: 'standard', 'wild'
spec.uncon    = 'true';         % options: 'true', 'false'
spec.paradigm = 'frequentist';  % options: 'Bayesian', 'frequentist'

for mm = 1:length(est)

    % LOAD MODEL SPECIFICATION
    est{mm}.model = spec.models(mm);     	% model
    est{mm}.id    = spec.id(mm);           % identification
    est{mm}.n_lag = spec.n_lag(mm);        % number of lags

    % VARIABLE TO NORMALIZE IRF TO
    % >> matters only for plotting IRF
    % >>  -1: 1 std dev; only for internal IV VAR;
    % >>   0: 1 unit change in instrument
    % >> i>0: 1 unit change in endogenous variable i
    est{mm}.norm_i = 0;

end


end