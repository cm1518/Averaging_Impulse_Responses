function pd_out = OLSfreq_pd_btstrp(Y, X, est_out, test_data, spec, n_lag, bb)

% ========================================================================
% COMPUTE PREDICTIVE DENSITY FROM OLS, varargin
% ========================================================================
%
% -- MODEL --
%   Y = X*beta + e
%
% -- INPUT --
%   Y (T-by-1): dependent variable
%   X (T-by-N): independent variables
%   b_hat, V: asymptotic distribution sqrt(T)*(beta-b_hat) ~ N(0, V*T)
%   train: observations to use for bootstrap errors
%   test: observations to compute predictive density for
%   n_draw: number of draws (same for parameter and shock draws)
%
% -- OUTPUT --
%   pd: predictive density for each observation in smpl
%
% ========================================================================


% ADD CONSTANT IF NEEDED
% PARSE INPUTS
% p = inputParser;
% addOptional(p, 'constant', 'true', ...
%     @(x) any(validatestring(x,{'true', 'false'})));
% addOptional(p, 'se', 'homoskedastic',...
%     @(x) any(validatestring(x,{'homoskedastic','white','nw'})));
% parse(p, varargin{:})




% ------------------------------------------------------------------------
% PRELIMINARIES
% ------------------------------------------------------------------------

% INCLUDE CONSTANT IF NEEDED
if strcmp(est_out.constant,'true')
    X = [X, ones(size(X,1),1)];
end

% TRAINING SAMPLE
Y_train = Y(est_out.sample);
X_train = X(est_out.sample,:);

% REMOVE NaN AG: I don't think this is necessary?
no_NaN = all(~isnan([Y_train,X_train]),2);
Y_train = Y_train(no_NaN,:);
X_train = X_train(no_NaN,:);

% TESTING SAMPLE (AG ADDED)
Y_test = test_data.Y;
if size(test_data.X,2) == size(X_train,2)
    X_test = test_data.X;
else
    X_test = test_data.X(:,2:end);
end

% LENGTH OF TRAINING AND TEST SAMPLE
T_train = length(Y_train);
T_test  = length(Y_test);

% DO BLOCK BOOTSTRAP FOR NEWEY-WEST STANDARD ERRORS
if strcmp(est_out.se,'nw')
    blk = floor(0.75*T_train^(1/3)); % Andrews (1991)
    % blk = floor(1.3*sqrt(T_train)); % Lazarus, Lewis, Stock, Watson (2018)
else
    blk = 1;
end


% ------------------------------------------------------------------------
% PREDICTIVE DENSITY
% ------------------------------------------------------------------------

% 0. LOAD PARAMETER ESTIMATES
b_hat = est_out.B;
V_hat = est_out.V;

% 1. DRAW PARAMETERS AND COMPUTE FITTED SHOCKS
b_draw = mvnrnd(b_hat, V_hat, spec.n_draw);
e_hat = Y_train - X_train*b_draw';

% 2. SAMPLE FITTED SHOCKS
if blk == 1
    e_hat = reshape(e_hat, spec.n_draw*T_train, 1);
    r = randi(T_train*spec.n_draw, T_test, spec.n_draw);
    e_rnd = e_hat(r);
else
    TT = (T_train-blk+1)*spec.n_draw;
    % matrix of blocks of shocks
    e_lag = NaN(blk, TT);
    for i = 1:blk
        e_lag(i,:) = reshape(e_hat(i:(T_train-blk)+i,:), 1, TT);
    end
    % randomize across blocks
    r_blk = randi(TT, ceil(T_test/blk)+1, spec.n_draw);
    e_blk = reshape(e_lag(:,r_blk), (ceil(T_test/blk)+1)*blk,  spec.n_draw);
    % randomize starting index
    r_st = randi(blk, 1, spec.n_draw);
    e_rnd = NaN(T_test, spec.n_draw);
    for i = 1:spec.n_draw
        e_rnd(:,i) = e_blk(r_st(i):r_st(i)+T_test-1, i);
    end
end

% 4. FOR EACH DRAW, COMPUTE DENSITY
xb = Y_test - e_rnd;
pd_draw = normpdf( xb, ...
    X_test*b_hat, sqrt(sum(X_test*V_hat.*X_test,2)) );
pd = mean(pd_draw,2)+eps^2;
if (1<bb) && (bb<(spec.n_blk+1))
    pd = pd(n_lag+1:end);
end
if bb<spec.n_blk
    pd = pd(1:end-(spec.H+1)); % AG added to make the structure build correctly
end
pd_out.pd = pd;


% ------------------------------------------------------------------------
% POINT FORECAST
% ------------------------------------------------------------------------

Y_hat = X_test*b_hat;
if (1<bb) && (bb<(spec.n_blk+1))
    Y_hat = Y_hat(n_lag+1:end);
end
if bb<spec.n_blk
    Y_hat = Y_hat(1:end-(spec.H+1)); % AG added to make the structure build correctly
end
pd_out.Y_hat = Y_hat;



end