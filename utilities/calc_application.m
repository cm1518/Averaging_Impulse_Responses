% ========================================================================
% EMPIRICAL APPLICATION: MAIN CALCULATIONS FOR AVERAGING IMPULSE RESPONSES
% ------------------------------------------------------------------------
% > estimate impulse responses
% > compute predictive densities
% > compute optimal weights
% > extract and save output of interest
% ========================================================================

fprintf('COMPUTATIONS STARTED: %s.\n', datestr(now))
delete(gcp('nocreate'))         % exit current parallel pool
parpool(feature('numcores'))	% new parallel pool with all available cores
tic;

% ESTIMATION
fprintf('estimation and predictive densities...\n')
parfor mm = 1:length(est)
    if strcmp(spec.paradigm,'Bayesian')
        if strcmp(spec.models{mm},'VAR')
            est{mm} = VARbayes_full(data, est{mm}, spec);
        elseif contains(spec.models{mm},'LP')
            est{mm} = LPbayes_full(data, est{mm}, spec);
        end
    elseif strcmp(spec.paradigm,'frequentist')
        if strcmp(spec.models{mm},'VAR')
            est{mm} = VARfreq_full(data, est{mm}, spec);
        elseif strcmp(spec.models{mm},'LP')
            est{mm} = LPfreq_full(data, est{mm}, spec);
        elseif strcmp(spec.models{mm},'FAVAR')
            est{mm} = FAVARfreq_full(data, est{mm}, spec);
        elseif contains(spec.models{mm},'single equation')
            est{mm} = sgleqn_full(data, est{mm}, spec);
        end
    end
end

% OPTIMAL WEIGHTS
fprintf('optimal weights...\n')
wt = wt_full(data, est, spec, var_list);    % compute weights
wt = avg_irf(est, spec, wt);                % average IRFs

% SAVE OUTPUT
fprintf('saving output...\n')
save(strcat('out/',out_name,'_',datestr(date,'yyyy-mm-dd')),'data','spec','est','wt')
fprintf('COMPLETED IN %.2f MINUTES.\n', toc/60)
