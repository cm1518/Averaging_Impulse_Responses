function est = LPbayes_full(data, est, spec)


% ========================================================================
% PRELIMINARIES
% ========================================================================

est.data.Y = table2array(data.tab(:,est.var));      % endogenous variables
[T,N] = size(est.data.Y);                           % no. of periods/endogenous variables

est.data.m = table2array(data.tab(:,est.IV));  	% instrument
est.Y_VAR  = [est.data.m, est.data.Y];          % VAR to be estimated
shock_i    = 1;                                 % shock of interest
endog_i    = 2:N+1;                             % endogenous variable indices
norm_i     = est.norm_i + 1;                    % IRF to normalize to


% ========================================================================
% ESTIMATION AND PREDICTIVE DENSITIES
% >> LOOP OVER BLOCKS
% ========================================================================

pd = cell(spec.H+1, 1);

for bb = 1:spec.n_blk+1

    % SUBSAMPLES
    test  = return_test_index(data.T, spec.n_blk, bb);
    Y_blk = est.Y_VAR;
    if bb <= spec.n_blk
        Y_blk(test,:) = NaN;
    end

    % ESTIMATION
    prior = VARbayes_prior_naive(Y_blk);
    out = LPbayes_est(Y_blk, est.n_lag, prior, spec.H, 'constant', spec.constant);
    
    % POSTERIOR DRAWS
    draw = LPbayes_draw(out, spec.n_draw);

    % CALCULATE PREDICTIVE DENSITIES
    pd = LPbayes_pd(est.data, draw, test, endog_i, spec.H, est.id);
%     est.pd{bb}.pd = NaN(length(test),N,spec.H);
%     for ii = 1:N
%         pd = LPbayes_pd(est.data, draw, test, ii, spec.H, est.id);
%         est.pd{bb}.pd(:,ii,:) = pd.pd;
%     end
    
    % COMPUTE IRFS
    irf = LPbayes_irf(draw, 1);                                            % 1 std dev shock
    
    % SAVE OUTPUT
    est.pd{bb} = pd;
    if norm_i >= 1
        est.irf{bb} = irf(end-N+1:end,:,:,:)./irf(norm_i,1,1,:);            % 1 unit shock
    else
        est.irf{bb} = irf(end-N+1:end,:,:,:);                               % 1 std dev shock
    end
    
end


end