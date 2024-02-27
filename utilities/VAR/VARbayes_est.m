function post = VARbayes_est(Y_data, n_lag, prior, varargin)

% ========================================================================
% BAYESIAN VAR ESTIMATION WITH NORMAL-INVERSE WISHART PRIOR
% ========================================================================
%
% -- MODEL --
%   y(t) = cons + B*[y(t-1),...,y(t-n_lag)] + e(t)
%
% -- INPUT --
%   Y(T-by-N): data
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
addOptional(p, 'sample', 1:size(Y_data,1), ...
    @(x) isnumeric(x) && min(floor(x)==x));
parse(p, varargin{:})
post = p.Results;


% ------------------------------------------------------------------------
% SET UP DATA MATRICES
% ------------------------------------------------------------------------

% SELECT SAMPLE
Y = Y_data(post.sample,:); [T,~] = size(Y);

% CONSTRUCT X WITH LAGGED Y
if strcmp(post.constant,'true')
    X_data = [lagmatrix(Y,1:n_lag),ones(T,1)];    % include constant in last column
else
    X_data = lagmatrix(Y,1:n_lag);                % omit constant
end

% REMOVE NaN
no_NaN = all(~isnan([Y_data,X_data]),2);
Y = Y_data(no_NaN,:); X = X_data(no_NaN,:);
[T,N] = size(Y);
post.sample_no_NaN = post.sample(no_NaN);


% ------------------------------------------------------------------------
% PRIOR HYPERPARAMTERS
% ------------------------------------------------------------------------

psi = prior.psi;    % prior mean of SIGMA
lam = prior.lam;    % MN prior: Inf => uninformative
alp = prior.alp; 	% MN prior: relative variance of distant lags
muu = prior.muu;    % sum-of-coefficients prior: Inf => uninformative, 0 => unit root + no cointegration
del = prior.del;	% dummy-initial-observation prior: Inf => uninformative


% ------------------------------------------------------------------------
% SET UP PRIOR
% ------------------------------------------------------------------------

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
        (d-N-1)*(lam^2)*(1/(s^alp))./psi_prior;	% covariance for lag coefficients
end
if strcmp(post.constant,'true')
    omega(end) = 1e3;	% variance for constant term
end

% 2. SUM-OF-COEFFICIENTS PRIOR
if muu < Inf
    Y_dum_1 = diag(mean(Y_data(1:n_lag,:),1))/muu;
    X_dum_1 = [repmat(Y_dum_1,1,n_lag),zeros(N,1)];
    Y = [Y_dum_1;Y]; X = [X_dum_1;X];
    T = T + N;
end

% 3. DUMMY-INITIAL-OBSERVATION PRIOR
if del < Inf
    Y_dum_2 = mean(Y_data(1:n_lag,:),1)/del;
    X_dum_2 = [repmat(Y_dum_2,1,n_lag),1/del];
    Y = [Y_dum_2;Y]; X = [X_dum_2;X];
    T = T + 1;
end


% ------------------------------------------------------------------------
% COMPUTE POSTERIOR
% ------------------------------------------------------------------------

% MEAN OF VAR COEFFICIENTS
B_hat = (X'*X + diag(1./omega))\(X'*Y + diag(1./omega)*B_prior);

% VAR RESIDUALS
e_hat = Y - X*B_hat;

% MODE OF COVARIANCE MATRIX
PSI_hat = ( e_hat'*e_hat + diag(psi) + ...
    (B_hat-B_prior)'*diag(1./omega)*(B_hat-B_prior) );

% LOG MARGINAL LIKELIHOOD
if muu == Inf && del == Inf
    
    % matrices to compute eigenvalues for
    % (more numerically stable than computing determinants)
    aaa = diag(sqrt(omega))*(X'*X)*diag(sqrt(omega));
    bbb = diag(1./sqrt(psi))*(e_hat'*e_hat + (B_hat-B_prior)'*diag(1./omega)*(B_hat-B_prior))*diag(1./sqrt(psi));
    
    % eigenvalues
    eigaaa = real(eig(aaa));	eigaaa(eigaaa<1e-12) = 0;
    eigbbb = real(eig(bbb));	eigbbb(eigbbb<1e-12) = 0;
    
    % log marginal likelihood
    logML = -N*T*log(pi)/2 + sum(gammaln((T+d-(0:N-1))/2)-gammaln((d-(0:N-1))/2)) ...
        - T*sum(log(psi))/2 - N*sum(log(eigaaa+1))/2 - (T+d)*sum(log(eigbbb+1))/2;
    
end


% ------------------------------------------------------------------------
% SAVE OUTPUT
% ------------------------------------------------------------------------

post.B_cons = B_hat;
post.cons = B_hat(end,:)';
post.B = B_hat(1:end-1,:)';
post.OMEGA = eye(N*n_lag+1)/(diag(1./omega) + X'*X);
post.PSI = PSI_hat;
post.df = T + d;
post.logML = logML;


end