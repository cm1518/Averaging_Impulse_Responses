function est = VARfreq_full(data, est, spec)

% ========================================================================
% PRELIMINARIES
% ========================================================================

est.data.Y = table2array(data.tab(:,est.var));      % endogenous variables
[T,N] = size(est.data.Y);                           % no. of periods/endogenous variables
if strcmp(est.id,'Cholesky')
    est.Y_VAR  = est.data.Y;                        % VAR to be estimated
    shock_i    = N;                                 % shock of interest
    endog_i    = 1:N;                               % endogenous variable indices
    norm_i     = est.norm_i;                        % IRF to normalize to
elseif strcmp(est.id,'internal IV')
    est.data.m = table2array(data.tab(:,est.IV));  	% instrument
    est.Y_VAR  = [est.data.m, est.data.Y];          % VAR to be estimated
    shock_i    = 1;                                 % shock of interest
    endog_i    = 2:N+1;                             % endogenous variable indices
    norm_i     = est.norm_i + 1;                    % IRF to normalize to
end


% ========================================================================
% ESTIMATION AND PREDICTIVE DENSITIES
% >> LOOP OVER BLOCKS
% ========================================================================

for bb = 1:spec.n_blk+1
    
    % SUBSAMPLES
    test  = return_test_index(data.T, spec.n_blk, bb);
    Y_blk = est.Y_VAR;
    if bb <= spec.n_blk
        Y_blk(test,:) = NaN;
    end

    % ESTIMATION
    out = VARfreq_est(Y_blk, est.n_lag, 'sample', 1:T, 'constant', spec.constant);
    
    % BOOTSTRAP DRAWS
    draw = VARfreq_btstrp(Y_blk, out, spec.n_draw, 'type', spec.btstrp);

    % CALCULATE PREDICTIVE DENSITIES
    % conditional on shock
    pd = VAR_pd(est.data, draw, test, endog_i, spec.H, est.id);
    % unconditional
    if strcmp(spec.uncon,'true')
        pd_uncon = VAR_pd_uncon(est.data, draw, test, endog_i, spec.H, est.id);
        est.pd_uncon{bb} = pd_uncon;
    end

    % COMPUTE IRFS
    irf = VAR_irf(draw, spec.H, shock_i, est.id, est.data);                 % 1 std dev shock
    
    % SAVE OUTPUT
    est.pd{bb} = pd;
    if norm_i >= 1
        est.irf{bb} = irf(end-N+1:end,:,:,:)./irf(norm_i,1,1,:);            % 1 unit shock
    else
        est.irf{bb} = irf(end-N+1:end,:,:,:);                               % 1 std dev shock
    end
    
end

end