function sim = sel_out_sim_SW(est, wt, spec)

for mm = 1:spec.n_model+1

    if mm <= spec.n_model
        irf_draw = squeeze(est{mm}.irf{end});
    elseif mm == spec.n_model + 1
        irf_draw = squeeze(wt.irf_av);
%     elseif mm == spec.n_model + 2
%         irf_draw = squeeze(wt_LS.irf_av);
    end

    % MEAN AND STANDARD DEVIATION
    sim.irf_mean{mm} = squeeze(mean(irf_draw,ndims(irf_draw),'omitnan'));
    sim.irf_std{mm}  = std(irf_draw,0,ndims(irf_draw),'omitnan');

end

% WEIGHTS
sim.wt = squeeze(wt.wt_blk);

end