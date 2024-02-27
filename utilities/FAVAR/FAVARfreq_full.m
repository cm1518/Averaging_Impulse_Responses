function est = FAVARfreq_full(data, est, spec)

% ========================================================================
% PRELIMINARIES
% ========================================================================

est.data.fac = table2array(data.tab(:,est.fac));    % endogenous variables
est.data.var = table2array(data.tab(:,est.var));    % endogenous variables
[T,N] = size(est.data.fac);                         % no. of periods/factors
est.Y = est.data.var;                               % variables of interest
n_var  = size(est.Y, 2);
if strcmp(est.id,'Cholesky')
    est.F_VAR  = est.data.fac;                      % VAR to be estimated
    shock_i    = N;                                 % shock of interest
elseif strcmp(est.id,'internal IV')
    est.data.m = table2array(data.tab(:,est.IV));  	% instrument
    est.F_VAR  = [est.data.m, est.data.fac];        % VAR to be estimated
    shock_i    = 1;                                 % shock of interest
end



% ========================================================================
% ESTIMATION AND PREDICTIVE DENSITIES
% >> LOOP OVER BLOCKS
% ========================================================================

for bb = 1:spec.n_blk+1
    
    % SUBSAMPLES
    test  = return_test_index(data.T, spec.n_blk, bb);
    F_blk = est.F_VAR;
    Y_blk = est.Y;
    if bb <= spec.n_blk
        F_blk(test,:) = NaN;
        Y_blk(test,:) = NaN;
        pd = NaN(length(test),n_var,spec.H+1);
    else
        pd = NaN(T,n_var,spec.H+1);
    end
    if strcmp(spec.uncon,'true')
        pd_uncon = pd;
    end
    

    % --------------------------------------------------------------------
    % VAR FOR FACTORS
    % --------------------------------------------------------------------

    % ESTIMATION
    out_F = VARfreq_est(F_blk, est.n_lag, 'sample', 1:T, 'constant', spec.constant);
    
    % BOOTSTRAP DRAWS
    draw_F = VARfreq_btstrp(F_blk, out_F, spec.n_draw, 'type', spec.btstrp);

    % DRAW FACTORS FROM PREDICTIVE DENSITY
    % conditional
    F_draw = FAVAR_pd_fac(est.data, draw_F, test, spec.H, est.id);
    % unconditional
    if strcmp(spec.uncon,'true')
        F_draw_uncon = FAVAR_pd_fac_uncon(est.data, draw_F, test, spec.H, est.id);
    end

    % IMPULSE RESPONSE
    irf_F = VAR_irf(draw_F, spec.H, shock_i, est.id, est.data);                 % 1 std dev shock

    
    % --------------------------------------------------------------------
    % LINEAR REGRESSION FOR VARIABLES OF INTEREST
    % --------------------------------------------------------------------
    
    out_Y  = cell(n_var, 1);
    irf_Y  = NaN( n_var, 1, spec.H+1, spec.n_draw);
    
    for ii = 1:size(Y_blk,2)

        % ESTIMATION
        out_Y{ii} = OLSfreq_est(Y_blk(:,ii), F_blk, 'se', 'white', 'constant', spec.constant);

        % PREDICTIVE DENSITIES
        % conditional
        pd(:,ii,:) = FAVAR_pd_var(Y_blk, est.Y(test,:), ...
            F_blk, F_draw, out_Y{ii}, spec);
        % unconditional
        if strcmp(spec.uncon,'true')
            pd_uncon(:,ii,:) = FAVAR_pd_var(Y_blk, est.Y(test,:), ...
                F_blk, F_draw_uncon, out_Y{ii}, spec);
        end
        
        % IMPULSE RESPONSE
        irf_Y(ii,1,:,:) = FAVAR_irf_draw(irf_F, out_Y{ii});
        
    end

    
    % --------------------------------------------------------------------
    % SAVE OUTPUT
    % --------------------------------------------------------------------
    
    est.pd{bb}.pd = pd;
    if strcmp(spec.uncon,'true')
        est.pd_uncon{bb}.pd = pd_uncon;
    end
    
    if est.norm_i > 0
        est.irf{bb} = irf_Y./irf_F(est.norm_i,1,1,:);   % 1 unit shock
    else
        est.irf{bb} = irf_Y;                            % 1 std dev shock
    end
    
    
end

end