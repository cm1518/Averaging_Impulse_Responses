% ========================================================================
% SELECT SPECIFICATIONS AND MODELS TO BE ESTIMATED
% >> SECTION 4.2: TFP SHOCK APPLICATION
% ========================================================================

% MODELS
% options: 'VAR', 'LP', 'FAVAR', 'single equation'
spec.models = {'VAR', 'VAR', 'VAR', 'FAVAR', ...
    'LP', 'LP', 'LP', 'single equation'};

% IDENTIFICATION SCHEMES
% >> VAR: 'Cholesky', 'internal IV'
% >> LP: 'recursive', 'non-recursive'
% >> single equation: N/A, leave as ''
spec.id = {'Cholesky', 'Cholesky', 'Cholesky', 'Cholesky', ...
    'recursive', 'recursive', 'recursive', ''};

% NUMBER OF LAGS, MODELS
spec.n_lag = 4*ones(length(spec.models),1);
spec.n_model = size(spec.models,2); 

% SELECT VARIABLES FOR EACH MODEL
est = cell(spec.n_model,1);
% >> 1. VAR: Cholesky, all controls
est{1}.var = {'unemp','GDP growth','utilization','fac_1','fac_2','fac_3','TFP growth'};
% >> 2. VAR: Cholesky, control for 1st 3 PCs
est{2}.var = {'unemp','fac_1','fac_2','fac_3','TFP growth'};
% >> 3. VAR: Cholesky, no controls
est{3}.var = {'unemp','TFP growth'};
% >> 5. FAVAR
est{4}.fac = {'fac_1', 'fac_2', 'fac_3', 'TFP growth'};
est{4}.var = {'unemp'};
% >> 5. LP: recursive, all controls
est{5}.var = {'unemp','GDP growth','utilization','fac_1','fac_2','fac_3'};
est{5}.IV  = 'TFP growth';
% >> 6. LP: recursive, control for 1st 3 PCs
est{6}.var = {'unemp','fac_1','fac_2','fac_3'};
est{6}.IV  = 'TFP growth';
% >> 7. LP: recursive, no controls
est{7}.var = {'unemp'};
est{7}.IV  = 'TFP growth';
% >> 8. single equation
est{8}.var = {'unemp'};
est{8}.IV  = 'TFP growth';


% ------------------------------------------------------------------------
% ESTIMATION DETAILS
% ------------------------------------------------------------------------

% ESTIMATION PARAMETERS
spec.H        = 20;             % number of horizons
spec.n_draw   = 2.5e3;          % number of draws
spec.n_blk    = 2;              % number of blocks for sample splitting. 0 indicates you train and test on the entire dataset
spec.constant = 'true';         % options: 'true', 'false'
spec.btstrp   = 'wild';         % options: 'standard', 'wild'
spec.uncon    = 'true';         % options: 'true', 'false'
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
    if mm <= 3
        est{mm}.norm_i = length(est{mm}.var);   % normalize IRF of last variable to be 1 on impact
    elseif mm == 4
        est{mm}.norm_i = length(est{mm}.fac);   % normalize to last factor
    else
        est{mm}.norm_i = 0;                     % LPs normalized automatically to instrument
    end
    
end

