function sim = sel_out_sim_HJ(est, wt, spec)

for mm = 1:spec.n_model+1

    if mm <= spec.n_model
        irf_draw = squeeze(est{mm}.irf{end}(1,:,:,:));
    elseif mm == spec.n_model + 1
        irf_draw = squeeze(wt.irf_av);
    end

    % MEAN AND STANDARD DEVIATION
    sim.irf_mean{mm} = squeeze(mean(irf_draw,2,'omitnan'));
%     sim.irf_std{mm}  = std(irf_draw,0,2,'omitnan');
    
    % QUANTILE THAT TRUE IRF IS IN
    sim.irf_qtl{mm} = sum(irf_draw<=sim_par.rho.^(0:20)',2)/size(irf_draw,2);
    
end

% WEIGHTS
sim.wt = squeeze(wt.wt_blk);

end