function est = LPfreq_full(data, est, spec)


% ========================================================================
% PRELIMINARIES
% ========================================================================

n_var = length(est.var);

% COLLECT DATA
est.data.Y = table2array(data.tab(:,est.var));
est.data.m = table2array(data.tab(:,est.IV));
est.Y_LP = [est.data.m, est.data.Y];


% ========================================================================
% ESTIMATION AND PREDICTIVE DENSITIES
% >> LOOP OVER BLOCKS
% ========================================================================

for bb = 1:spec.n_blk+1
    
    % --------------------------------------------------------------------
    % TEST AND TRAINING SUBSAMPLES
    % --------------------------------------------------------------------
    
    test  = return_test_index(data.T, spec.n_blk, bb);

    if bb > 1 && bb < spec.n_blk+1
        start_index = test(1)-est.n_lag;
    else
        start_index = test(1);
    end

    if bb < spec.n_blk
        Y_dropped = est.data.Y(start_index:test(end)+spec.H+1,:);
        z_dropped = est.data.m(start_index:test(end)+spec.H+1,:);
    elseif bb >= spec.n_blk
        Y_dropped = est.data.Y(start_index:test(end),:);
        z_dropped = est.data.m(start_index:test(end),:);
    end

    % Make Test Block NaN (for non global case)
    Y_tr_smpl = est.data.Y;
    z_tr_smpl = est.data.m;
    if bb <= spec.n_blk && spec.n_blk~=1
        Y_tr_smpl(test,:) = NaN;
        z_tr_smpl(test,:) = NaN; 
    end

    
    % --------------------------------------------------------------------
    % COMPUTATIONS
    % --------------------------------------------------------------------
    
    out      = cell(spec.H+1, 1);
    pd       = cell(spec.H+1, 1);
    pd_uncon = cell(spec.H+1, 1);
    irf_tmp  = NaN(    n_var, 1, spec.H+1, spec.n_draw);
    for ii = 1:n_var        
        
        for h = 1:spec.H+1

            % LHS: VARIABLE OF INTEREST hor PERIODS AHEAD
            Y_h = lagmatrix(Y_tr_smpl(:,ii), -h+1);

            % RHS
            if strcmp(est.id,'non-recursive')
                % lagged controls only
                X = [ z_tr_smpl, lagmatrix([z_tr_smpl,Y_tr_smpl], 1:est.n_lag)];    
            elseif strcmp(est.id,'recursive')
                % include contemporaneous controls
                endog_var_i     = 1:n_var;
                endog_var_i(ii) = [];       % don't include current LHS var
                X = [ z_tr_smpl, Y_tr_smpl(:,endog_var_i), lagmatrix([z_tr_smpl,Y_tr_smpl], 1:est.n_lag) ];
            elseif strcmp(est.id,'lag_endog')
                % lagged endogenous y, for Herbst-Johannsen AR(1) Monte Carlo
                X = [ z_tr_smpl, lagmatrix(Y_tr_smpl, 1:est.n_lag) ];
            end
            
            % ESTIMATE LP
            out{h} = OLSfreq_est(Y_h, X, 'se', 'white', 'constant', spec.constant);
            est.out{bb}.B(:,ii,h) = out{h}.B;
            est.out{bb}.V(:,:,ii,h) = out{h}.V;

            % PREDICTIVE DENSITIES
            % test sample: only used in LP Freq because we cannot use data_est
            test_data = create_test_data(spec, est, Y_dropped, z_dropped, ii, h);
            % predictive density
            pd{h} = OLSfreq_pd_btstrp(Y_h, X, out{h}, test_data, spec, est.n_lag, bb);
            est.pd{bb}.pd(:,ii,h) = pd{h}.pd;
            est.pd{bb}.pt(:,ii,h) = pd{h}.Y_hat;
            
            % UNCONDITIONAL PREDICTIVE DENSITIES
            if strcmp(spec.uncon,'true')
                out_uncon = OLSfreq_est(Y_h, X(:,2:end), 'se', 'white', 'constant', spec.constant);
                pd_uncon{h} = OLSfreq_pd_btstrp(Y_h, X(:,2:end), out_uncon, test_data, spec, est.n_lag, bb);
                est.pd_uncon{bb}.pd(:,ii,h) = pd_uncon{h}.pd;
            end

        end

        % COMPUTE IRFS
        irf_tmp(ii,1,:,:)         = LPfreq_irf_draw(out, spec.n_draw);
        est.irf_qtl{bb}(ii,1,:,:) = LPfreq_irf(out, normcdf([-1,0,1]));

    end
    
    % NORMALIZE IRFS TO UNIT RESPONSE ON IMPACT OF SELECTED VARIABLE
    if est.norm_i > 0
        est.irf{bb} = irf_tmp./irf_tmp(est.norm_i,1,1,:);
    else
        est.irf{bb} = irf_tmp;
    end
    
end


end