function wt = wt_full_LS(data, est, spec, var)

% VARIABLES TO COMPUTE WEIGHTS FOR
wt.var = var;                     

% COMPUTE OPTIMAL WEIGHTS
% >> conditional: always compute
wt = wt_subfn(data, est, spec, wt); %, 'false');
% % >> unconditional: if needed
% if strcmp(spec.uncon,'true')
%     wt_uncon = wt_subfn(data, est, spec, wt, 'true');
%     wt.wt_uncon_blk = wt_uncon.wt_blk;
%     wt.wt_uncon_glb = wt_uncon.wt_glb;
% end

% AVERAGE IRFS ACROSS MODELS
wt = avg_irf(est, spec, wt);

end

function wt = wt_subfn(data, est, spec, wt) %, uncon)

% ========================================================================
% CONTAINERS
% ------------------------------------------------------------------------
% blk: with sample split into blocks
% glb: without splitting sample into blocks
% ========================================================================

pt_stack_blk = cell(spec.n_model, length(wt.var));
pt_stack_glb = cell(spec.n_model, length(wt.var));
wt.wt_blk    =  NaN(spec.n_model, length(wt.var), spec.H+1);
wt.wt_glb    =  NaN(spec.n_model, length(wt.var), spec.H+1);


% ========================================================================
% OPTIMAL WEIGHTS
% ========================================================================

for ii = 1:length(wt.var)
    
    % COMBINE PREDICTIVE DENSITIES ACROSS BLOCKS
    for mm = 1:spec.n_model
        vv = find(strcmp(est{mm}.var, wt.var(ii)));
        % benchmark: with sample split into blocks
        pt_stack_blk{mm,ii} = NaN(data.T,spec.H+1);
        for bb = 1:spec.n_blk
            test = return_test_index(data.T, spec.n_blk, bb);
%             if strcmp(uncon,'true')
%                 pd_stack_blk{mm,ii}(test,:) = est{mm}.pd_uncon{bb}.pd(:,vv,:);
%             else
                pt_stack_blk{mm,ii}(test,:) = est{mm}.pd{bb}.pt(:,vv,:);
%             end
        end
        % global: without splitting sample into blocks
%         if strcmp(uncon,'true')
%             pd_stack_glb{mm,ii} = est{mm}.pd_uncon{spec.n_blk+1}.pd(:,vv,:);
%         else 
            pt_stack_glb{mm,ii} = est{mm}.pd{spec.n_blk+1}.pt(:,vv,:);
%         end
    end

    % CALCULATE WEIGHTS
    Y = table2array(data.tab(:,wt.var(ii)));
    for hh = 1:spec.H+1
        pt_blk_h = cell(spec.n_model,1);
        pt_glb_h = cell(spec.n_model,1);
        for mm = 1:spec.n_model
            pt_blk_h{mm} = pt_stack_blk{mm,ii}(:,hh);
            pt_glb_h{mm} = pt_stack_glb{mm,ii}(:,hh);
        end
        wt.wt_blk(:,ii,hh) = wt_opt_LS(Y, pt_blk_h);
        wt.wt_glb(:,ii,hh) = wt_opt_LS(Y, pt_glb_h);
    end
    
end



end