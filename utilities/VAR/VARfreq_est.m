function out = VARfreq_est(Y, n_lag, varargin)

% ========================================================================
% FREQUENTIST ESTIMATION OF VECTOR AUTOREGRESSION
% ========================================================================
%
% -- MODEL --
%   y(t) = con + B*[y(t-1),...,y(t-nlag)] + e
%
% -- INPUT --
%   Y (T-by-N): data
%   nlag: number of lags in X
%   inc_cons (optional): 0 to omit constant; include constant by default
%
% -- OUTPUT --
%   cons: estimate of constant term
%   B_hat: estimate of coefficient
%   V: covariance matrix for residuals
%   e: fitted residuals
%
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% PARSE INPUTS
% ------------------------------------------------------------------------

p = inputParser;
addOptional(p, 'constant', 'true', ... 
    @(x) any(validatestring(x,{'true', 'false'})));
addOptional(p, 'sample', 1:size(Y,1), ...
    @(x) isnumeric(x) && min(floor(x)==x));
parse(p, varargin{:})
out = p.Results;


% ------------------------------------------------------------------------
% SET UP DATA MATRICES
% ------------------------------------------------------------------------

% SELECT SAMPLE
Y = Y(out.sample,:); [T,N] = size(Y);

% CONSTRUCT X WITH LAGGED Y
if strcmp(out.constant,'true')
    X = [lagmatrix(Y,1:n_lag),ones(T,1)];    % include constant in last column
else
    X = lagmatrix(Y,1:n_lag);                % omit constant
end

% REMOVE ROWS THAT CONTAIN NaN
no_NaN = all(~isnan([Y,X]),2);
Y = Y(no_NaN,:); X = X(no_NaN,:);
out.sample_no_NaN = out.sample(no_NaN);


% ------------------------------------------------------------------------
% ESTIMATION
% ------------------------------------------------------------------------

B = X\Y;                                % coefficients
e = Y - X*B;                            % residuals
V = (e'*e)/(size(X,1) - size(X,2));     % variance
if strcmp(out.constant,'true')
    % include constant
    cons = B(end,:)';   B = B(1:end-1,:)';
else
    % omit constant
    cons = zeros(N,1);  B = B';
end


% ------------------------------------------------------------------------
% SAVE OUTPUT
% ------------------------------------------------------------------------

out.cons = cons;
out.B = B;
out.V = V;
out.e = e;

end