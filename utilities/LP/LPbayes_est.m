function post = LPbayes_est(data, n_lag, prior, H, varargin)

% ========================================================================
% BAYESIAN LINEAR PROJECTION ESTIMATION WITH NORMAL-INVERSE WISHART PRIOR
%   FOLLOWS MIRANDA-AGRIPPINO AND RICOO (2021)
% ========================================================================
%
% -- MODEL --
%   y(t) = cons + B*[y(t-1),...,y(t-n_lag)] + e(t)
%
% -- INPUT --
%   data:
%    >> Y(T-by-N): data
%    >> m(T-by-1): instrument
%   n_lag: number of lags
%   prior: conjugate prior hyperparameters
%   constant (optional): whether to include constant, 'true' by default
%   sample (optional): sample to use for estimation, full sample by default
%   
% -- OUTPUT --
%   post: parameters characterizing posterior
%
% ------------------------------------------------------------------------

% ------------------------------------------------------------------------
% PARSE INPUTS
% ------------------------------------------------------------------------

p = inputParser;
addOptional(p, 'constant', 'true', ... 
    @(x) any(validatestring(x,{'true', 'false'})));
%addOptional(p, 'sample', 1:size(data.Y,1), ...
%    @(x) isnumeric(x) && min(floor(x)==x));
parse(p, varargin{:})



post = cell(H+1,1);

for h = 0:H
    
    post{h+1} = p.Results;
    
    % --------------------------------------------------------------------
    % SET UP DATA MATRICES
    % --------------------------------------------------------------------
    
    % SELECT SAMPLE
    %Y_data = [data.m, data.Y];
    %Y = Y_data(post{h+1}.sample,:); 
    Y_data=data; Y=Y_data;
    [T,~] = size(Y);

    % CONSTRUCT X WITH LAGGED Y, INCLUDE INSTRUMENT IN CURRENT PERIOD
    if strcmp(post{h+1}.constant,'true')
        X_data = [lagmatrix(Y,1:n_lag),ones(T,1)];	% include constant in last column
    else
        X_data = lagmatrix(Y,1:n_lag);              % omit constant
    end
    
    % CONSTRUCT Y, h PERIODS AHEAD
    Y_data = lagmatrix(Y_data, -h);

    % REMOVE NaN
    no_NaN = all(~isnan([Y_data,X_data]),2);
    Y = Y_data(no_NaN,:); X = X_data(no_NaN,:);
    [T,N] = size(Y);
    %post{h+1}.sample_no_NaN = post{h+1}.sample(no_NaN);


    % --------------------------------------------------------------------
    % PRIOR HYPERPARAMTERS
    % --------------------------------------------------------------------

    psi = prior.psi;    % prior mean of SIGMA
    lam = prior.lam;    % MN prior: tightness
    alp = prior.alp; 	% MN prior: relative variance of distant lags
    muu = prior.muu;    % sum-of-coefficients prior: Inf => uninformative, 0 => unit root + no cointegration
    del = prior.del;	% dummy-initial-observation prior: Inf => uninformative


    % --------------------------------------------------------------------
    % SET UP PRIOR
    % --------------------------------------------------------------------

    % 0. INVERSE-WISHART FOR SIGMA
    d = N+2;          % degrees of freedom
    psi_prior = psi;  % diagonal of prior mean for SIGMA multiplied by (d-N-1)

    % 1. MINNESOTA PRIOR
    % >> mean
    B_prior = zeros(N*n_lag+1, N);
    B_prior(2:N,2:N) = eye(N-1);
    % >> variance
    omega = zeros(N*n_lag+1,1);
    for s = 1:n_lag
        omega((s-1)*N+1:s*N) = ...
            (d-N-1)*(lam^2)*(1/((s+h)^alp))./psi_prior;	% covariance for lag coefficients
    end
    if strcmp(post{h+1}.constant,'true')
        omega(end) = 1e3;	% variance for constant term
    end

    % 2. SUM-OF-COEFFICIENTS PRIOR
    if muu < Inf
        Y_dum_1 = diag(mean(Y(1:n_lag,:),1))/muu; %AG changed Y_data to Y to correct for NaN sample splitting
        X_dum_1 = [repmat(Y_dum_1,1,n_lag),zeros(N,1)];
        Y = [Y_dum_1;Y]; X = [X_dum_1;X];
        T = T + N;
    end

    % 3. DUMMY-INITIAL-OBSERVATION PRIOR
    if del < Inf
        Y_dum_2 = mean(Y(1:n_lag,:),1)/del; %AG changed Y_data to Y to correct for NaN sample splitting
        X_dum_2 = [repmat(Y_dum_2,1,n_lag),1/del];
        Y = [Y_dum_2;Y]; X = [X_dum_2;X];
        T = T + 1;
    end


    % --------------------------------------------------------------------
    % COMPUTE POSTERIOR
    % --------------------------------------------------------------------

    % MEAN OF VAR COEFFICIENTS
    B_hat = (X'*X + diag(1./omega))\(X'*Y + diag(1./omega)*B_prior);

    % VAR RESIDUALS
    e_hat = Y - X*B_hat;

    % MODE OF HAC COVARIANCE MATRIX
    HAC_lag = h-1;                                      % number of lags: set to 0 for non-HAC variance
    HAC_var = e_hat'*e_hat;                             % non-HAC variance (used for VAR)
    HAC_wts = (HAC_lag + 1 - (1:HAC_lag))./(HAC_lag+1); % weights for each lag in HAC variance
    for l = 1:HAC_lag
        GAM = e_hat(l+1:T,:)'*e_hat(1:T-l,:);
        HAC_var = HAC_var + HAC_wts(l)*(GAM+GAM');
    end
    PSI_hat = ( HAC_var + diag(psi) + ...
        (B_hat-B_prior)'*diag(1./omega)*(B_hat-B_prior) );



    % --------------------------------------------------------------------
    % SAVE OUTPUT
    % --------------------------------------------------------------------
    
    post{h+1}.B_cons = B_hat;
    post{h+1}.cons = B_hat(end,:)';
    post{h+1}.B = B_hat(1:end-1,:)';
    post{h+1}.OMEGA = eye(N*n_lag+1)/(diag(1./omega) + X'*X);
    post{h+1}.PSI = PSI_hat;
    post{h+1}.df = T + d;
    
end


end