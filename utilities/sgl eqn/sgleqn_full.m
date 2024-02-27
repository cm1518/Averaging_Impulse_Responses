function est = sgleqn_full(data, est, spec)


% ========================================================================
% PRELIMINARIES
% ========================================================================

% n_var = length(est.var);

% COLLECT DATA
est.data.Y = table2array(data.tab(:,est.var));
est.data.m = table2array(data.tab(:,est.IV));
est.dY_sgleqn = [NaN; diff(est.data.Y)];
est.X_sgleqn  = [lagmatrix(est.data.m,0:spec.H),...
    lagmatrix(est.dY_sgleqn,1:est.n_lag)];

% ========================================================================
% ESTIMATION AND PREDICTIVE DENSITIES
% >> LOOP OVER BLOCKS
% ========================================================================

for bb = 1:spec.n_blk+1
    
    % --------------------------------------------------------------------
    % TEST AND TRAINING SUBSAMPLES
    % --------------------------------------------------------------------
    
    % SUBSAMPLES
    test  = return_test_index(data.T, spec.n_blk, bb);
    dY_blk = est.dY_sgleqn;
    X_blk  = est.X_sgleqn;
    if bb <= spec.n_blk
        dY_blk(test,:) = NaN;
        X_blk(test,:)  = NaN;
    end

    
    % --------------------------------------------------------------------
    % COMPUTATIONS
    % --------------------------------------------------------------------

    % ESTIMATION
    out = OLSfreq_est(dY_blk, X_blk, 'constant', spec.constant);

    % BOOTSTRAP DRAWS
    btstrp = sgleqn_btstrp(dY_blk, X_blk, out, spec.n_draw, 'type', 'standard'); %spec.btstrp);

    % CALCULATE PREDICTIVE DENSITIES
    pd = sgleqn_pd(est, btstrp, test, spec.H);
    
    % UNCONDITIONAL PREDICTIVE DENSITIES
    if strcmp(spec.uncon,'true')
        out_uncon = OLSfreq_est(dY_blk, X_blk(:,2:end), 'constant', spec.constant);
        btstrp_uncon = sgleqn_btstrp(dY_blk, X_blk(:,2:end), out_uncon, spec.n_draw, 'type', 'standard');
        pd_uncon = sgleqn_pd(est, btstrp_uncon, test, spec.H);
        est.pd_uncon{bb} = pd_uncon;
    end
    
    % CONSTRUCT IRFs
    irf = sgleqn_irf_draw(btstrp, est.n_lag, spec.H);

    % SAVE OUTPUT
    est.pd{bb}  = pd;
    est.irf{bb} = NaN(1,1,size(irf,1),size(irf,2));
    est.irf{bb}(1,1,:,:) = repmat(irf,1,1);

end


end