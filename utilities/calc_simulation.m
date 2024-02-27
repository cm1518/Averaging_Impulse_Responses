% ========================================================================
% MONTE CARLO: MAIN CALCULATIONS FOR AVERAGING IMPULSE RESPONSES
% ------------------------------------------------------------------------
% > estimate impulse responses
% > compute predictive densities
% > compute optimal weights
% > extract and save output of interest
% ========================================================================

fprintf('MONTE CARLO STARTED: %s.\n', datestr(now))
delete(gcp('nocreate'))         % exit current parallel pool
parpool(feature('numcores'))	% new parallel pool with all available cores
tic;

% CELL TO STORE SIMULATION OUTPUT
sim = cell(n_sim,1);

% COMPUTATIONS
parfor ii = 1:n_sim
    
    rng(ii)
    
    % SIMULATE DATA
    data = sim_fn();

    % MODEL AND ESTIMATION SPECIFICATIONS
    [spec, est] = spec_fn();

    % ESTIMATES AND PREDICTIVE DENSITIES
    for mm = 1:length(est)
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
    wt = wt_full(data, est, spec, var_list);
    
    % STORE OUTPUT
    sim{ii} = sel_fn(est, wt, spec);

end

% PERMANENT MODEL SPECIFICATION STRUCTURES TO BE SAVED
[spec, est] = spec_fn();

% SAVE OUTPUT
save(strcat('out/',out_name,'_lean_',datestr(date,'yyyy-mm-dd')),...
    'sim_par','spec','est','sim','-v7.3')

% REPORT
fprintf('COMPLETED IN %.2f MINUTES.\n', toc/60)
