function prior = VARbayes_prior_naive(Y)

% ========================================================================
% SET DEFAULT PRIOR FOR VAR
% ========================================================================

% make prior mean of SIGMA same scale as data
prior.psi = 0.5*var(Y, 'omitnan'); % AG added omitnan to account for sample splitting	

% MN prior: tightness and relative variance of distant lags
prior.lam = 1e3;
prior.alp = 2;

% sum-of-coefficients prior: \infty => uninformative, 0 => unit root + no cointegration
prior.muu = Inf;

% dummy-initial-observation prior: \infty => uninformative
prior.del = Inf;

end