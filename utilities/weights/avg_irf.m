function wt = avg_irf(est, spec, wt)

if strcmp(spec.uncon,'true')
    wt.irf_av_uncon     = avg_irf_subfn(est, spec, wt, 'true');
    wt.irf_av_uncon_var = var(wt.irf_av_uncon,0,3);
end
wt.irf_av     = avg_irf_subfn(est, spec, wt, 'false');
wt.irf_av_var = var(wt.irf_av,0,3);

end

function irf_av = avg_irf_subfn(est, spec, wt, uncon)

% CONTAINER
irf_av = NaN(length(wt.var), spec.H+1, spec.n_draw);

% AVERAGE ACROSS MODEL
% sample from each model according to weights
for ii = 1:length(wt.var)
    
    % MAKE DRAWS
    for hh = 1:spec.H+1
        % draw across models based on optimal weights
        if strcmp(uncon,'true')
            smpl_a = mnrnd(spec.n_draw, wt.wt_uncon_blk(:,ii,hh));
        else
            smpl_a = mnrnd(spec.n_draw, wt.wt_blk(:,ii,hh));
        end
        % draw within models
        for mm = 1:spec.n_model
            i_st = find(isnan(squeeze(irf_av(ii,hh,:))),1);
            i_en = sum(smpl_a(1:mm));
            if i_en >= i_st
                vv = strcmp(est{mm}.var, wt.var{ii});
                irf_av(ii,hh,i_st:i_en) = datasample(...
                    est{mm}.irf{end}(vv,1,hh,:), smpl_a(mm));
            end
        end
    end
    
end

end