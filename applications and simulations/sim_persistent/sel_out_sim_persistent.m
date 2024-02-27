function sim = sel_out_sim_ARMA(est, wt, spec, sim_par)

for mm = 1:3

    if mm <= spec.n_model
        irf_draw = squeeze(est{mm}.irf{end});
    else
        irf_draw       = squeeze(wt.irf_av);
        irf_draw_uncon = squeeze(wt.irf_av_uncon);
    end

    % MEAN AND STANDARD DEVIATION
    sim.irf_mean{mm} = squeeze(mean(irf_draw,2,'omitnan'));
    sim.irf_std{mm}  = std(irf_draw,0,2,'omitnan');
    if mm > spec.n_model
        sim.irf_mean_uncon = mean(irf_draw_uncon,2,'omitnan');
        sim.irf_std_uncon  = std(irf_draw_uncon,0,2,'omitnan');
    end

    % QUANTILE THAT TRUE IRF IS IN
    sim.irf_qtl{mm} = sum(irf_draw<=sim_par.rho_y.^(0:20)',2)/size(irf_draw,2);

    % POWER
    sim.irf_pow{mm} = mean(irf_draw>0,2);

end

% WEIGHTS
sim.wt       = squeeze(wt.wt_blk);
sim.wt_uncon = squeeze(wt.wt_uncon_blk);

end