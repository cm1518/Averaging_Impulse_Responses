function btstrp = sgleqn_btstrp(dY, X, est_out, n_draw, varargin)

% ========================================================================
% BOOSTRAP DRAWS FOR VAR
% ========================================================================
%
% -- MODEL --
%   dy(t)     = cons + B*X + u(t)
%   var(u(t)) = V
%   X(t)      = [z(t),...,z(t-n_lag_z),dy(t-1),...,dy(t-n_lag)]
%
% -- INPUT --
%   dY(T-by-1), X(T-by-(n_lag_z+n_lag+1)): data
%   est_out: estimation output from OLSfreq_est.m
%   n_draw: number of bootstrap draws to make
%   type (optional): standard or wild; standard by default
%   
% -- OUTPUT --
% btstrp: structure with boostrap draws
%
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% PARSE INPUTS
% ------------------------------------------------------------------------

p = inputParser;
addOptional(p, 'type', 'standard', ... 
    @(x) any(validatestring(x,{'standard', 'wild'})));
parse(p, varargin{:})
type = p.Results.type;


% ------------------------------------------------------------------------
% PREPARE FOR BOOTSTRAP
% ------------------------------------------------------------------------

% LOAD OUTPUT FROM ESTIMATION
B_hat    = est_out.B;

% DIMENSIONS/SUBSAMPLES
dY = dY(est_out.sample,:);
dY_nonan = dY(est_out.sample_no_NaN,:);  % Added by AG so that the initial condition isn't NaN
X_nonan  =  X(est_out.sample_no_NaN,:);
T = length(dY_nonan);                  % Added by AG so that the loop runs for the correct size
if strcmp(est_out.constant,'true')
    X_nonan = [X_nonan, ones(T,1)];
end

% FITTED RESIDUALS
e_hat = dY_nonan - X_nonan*B_hat;



% ------------------------------------------------------------------------
% BOOTSTRAP
% ------------------------------------------------------------------------

% 1. DRAW RESIDUALS
if strcmp(type,'standard')
    r = randi(size(e_hat,1), T, n_draw);
elseif strcmp(type,'wild')
    r = randi([0,1], T, n_draw)*2 - 1;
end

% 2. COMPUTE BOOSTRAP TIME SERIES
if strcmp(type,'standard')
    e_btstrp = e_hat(r);
elseif strcmp(type,'wild')
    e_btstrp = e_hat.*r;
end
dY_btstrp = X_nonan*B_hat + e_btstrp;
btstrp.V = sum(e_btstrp.*e_btstrp)/(size(X_nonan,1)-size(X_nonan,2));

% 3. REESTIMATE PARAMETERS AND STORE DRAWS
btstrp.B = NaN(size(B_hat,1), n_draw);
for i = 1:n_draw
    out_tmp = OLSfreq_est(dY_btstrp(:,i), X_nonan, 'se', 'white', 'constant',  'false');
    btstrp.B(:,i)  = out_tmp.B;
end
btstrp.est_out = est_out;

end