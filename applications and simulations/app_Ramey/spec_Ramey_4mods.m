% ========================================================================
% SELECT SPECIFICATIONS AND MODELS TO BE ESTIMATED
% >> APPENDIX D.2: MONETARY SHOCK APPLICATION (WITH CHOLESKY VAR)
% ========================================================================

% MODELS
% >> 'VAR', 'LP', 'single equation', 'FAVAR'
spec.models = {'VAR', 'VAR', 'LP', 'LP'};

% IDENTIFICATION SCHEMES
% >> VAR: 'Cholesky', 'internal IV'
% >> LP: 'recursive', 'non-recursive'
% >> single equation: N/A, leave as ''
spec.id = {'Cholesky', 'internal IV', 'non-recursive', 'recursive'};

% NUMBER OF LAGS, MODELS
spec.n_lag = [12, 12, 2, 2];
spec.n_model = size(spec.models,2); 

% SELECT VARIABLES FOR EACH MODEL
est = cell(spec.n_model,1);
% >> 1. VAR: Cholesky
est{1}.var = {'LIP', 'UNEMP', 'LCPI', 'LPCOM', strcat('CUM',shock)};
% >> 2. VAR: internal instrument
est{2}.var = {'LIP', 'UNEMP', 'LCPI', 'LPCOM', 'FFR'};
est{2}.IV  = shock;
% >> 3. LP: non-recursive
est{3}.var = est{2}.var;
est{3}.IV  = est{2}.IV;
% >> 4. LP: recursive
est{4}.var = est{2}.var;
est{4}.IV  = est{2}.IV;


% ------------------------------------------------------------------------
% ESTIMATION DETAILS
% ------------------------------------------------------------------------

% ESTIMATION PARAMETERS
spec.H        = 48;             % number of horizons
spec.n_draw   = 5e3;            % number of draws
spec.n_blk    = 2;              % number of blocks for sample splitting. 0 indicates you train and test on the entire dataset 
spec.constant = 'true';         % options: 'true', 'false'
spec.btstrp   = 'wild';         % options: 'standard', 'wild'
spec.uncon    = 'false';        % options: 'true', 'false'
spec.paradigm = 'frequentist';  % options: 'Bayesian', 'frequentist'

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
    est{mm}.norm_i = length(est{mm}.var);   % normalize IRF of last variable (FFR) to be 1 on impact
    
end