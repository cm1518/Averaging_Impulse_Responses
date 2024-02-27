function btstrp = VARfreq_btstrp(Y, est_out, n_draw, varargin)

% ========================================================================
% BOOSTRAP DRAWS FOR VAR
% ========================================================================
%
% -- MODEL --
%   y(t) = cons + B*[y(t-1),...,y(t-n_lag)] + e(t)
%   var(u(t)) = V
%
% -- INPUT --
%   Y(T-by-N): data
%   est_out: estimation output from VARfreq_est.m
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
cons_hat = est_out.cons;
B_hat    = est_out.B;
e_hat    = est_out.e;

% DIMENSIONS
Y = Y(est_out.sample,:);
% [T,N] = size(Y); commented out by AG to sample split
Y_nonan = Y(est_out.sample_no_NaN,:); %Added by AG so that the initial condition isn't NaN
[T,N] = size(Y_nonan); %Added by AG so that the loop runs for the correct size
n_lag = size(B_hat,2)/N;

% MATRICES
btstrp.cons = NaN(N, n_draw);
btstrp.B    = NaN(N, N*n_lag, n_draw);
btstrp.V    = NaN(N, N, n_draw);


% ------------------------------------------------------------------------
% BOOTSTRAP
% ------------------------------------------------------------------------

% 1. DRAW RESIDUALS
if strcmp(type,'standard')
    r = randi(size(e_hat,1), T, n_draw);
elseif strcmp(type,'wild')
    r = randi([0,1], T-n_lag, n_draw)*2 - 1;
end

% 2. COMPUTE BOOSTRAP TIME SERIES
Y_btstrp = NaN(T, N*n_draw); 
Y_btstrp(1:n_lag,:) = repmat(Y_nonan(1:n_lag,:),1,n_draw); % use data for initial condition -- use nonan (AG)
for t = n_lag+1:T
    % shocks
    if strcmp(type,'standard')
        e_tmp = e_hat(r(t-n_lag,:),:)';
    elseif strcmp(type,'wild')
        e_tmp = repmat(e_hat(t-n_lag,:)',1,n_draw).*repmat(r(t-n_lag,:),N,1);
    end
    
    % lagged Y
    Y_lag = Y_btstrp(t-1:-1:t-n_lag,:);
    Y_lag = permute(reshape(Y_lag, n_lag, N, n_draw), [2,1,3]);
    Y_lag = reshape(Y_lag, n_lag*N, n_draw);
    
    % iterate forward
    Y_btstrp(t,:) = repmat(cons_hat, n_draw, 1)...
        + reshape(B_hat*Y_lag, N*n_draw, 1) + reshape(e_tmp, N*n_draw, 1);
end
Y_btstrp = reshape(Y_btstrp, T, N, n_draw);

% 3. REESTIMATE PARAMETERS AND STORE DRAWS
for i = 1:n_draw
    out_tmp = VARfreq_est(Y_btstrp(:,:,i), n_lag);
    btstrp.cons(:,i) = out_tmp.cons;
    btstrp.B(:,:,i)  = out_tmp.B;
    btstrp.V(:,:,i)  = out_tmp.V;
end
btstrp.est_out = est_out;

end