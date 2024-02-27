% ========================================================================
% PLOT FOR MONTE CARLO WITH PERSISTENT SHOCK
% ========================================================================

clear
clc
addpath(genpath(pwd))

% LOAD OUTPUT
load('out/sim_persistent_lean_2023-09-06.mat')

% PLOT COLORS
col = [76, 114, 176;
    221, 132, 82;
    85, 168, 104;
    129, 114, 179;
    196, 78, 82;
    204, 185, 116;
    218, 139, 195;
    100,181,205;
    0,120,140;
    170, 10, 0]/255;

% PLOT FIGURE 2
close all
ff = figure('position',[100,100,900,500]);
sl = 0;
lw = 1.5;
for ii = 1:4
    colororder(col)
    subplot('position', [mod(ii+1,2)*0.39+0.06, 1.05-ceil(ii/2)*0.47, 0.31, 0.33])
    if ii == 1
        wt_LP = mean(cell2mat(cellfun(@(x) x.wt(2,:)', sim, 'UniformOutput', false)),2);
        wt_LP_uncon = mean(cell2mat(cellfun(@(x) x.wt_uncon(2,:)', sim, 'UniformOutput', false)),2);
        plot(0:spec.H,wt_LP,'--','linewidth',lw,'color',[0,0,0])
        hold on
        plot(0:spec.H,wt_LP_uncon,'--','linewidth',lw,'color',col(end-1,:))
        ylim([0,Inf])
        title('Weight on LP')
    elseif ii == 2
        tmp = cell2mat(cellfun(@(x) x.wt(2,:)', sim, 'UniformOutput', false));
        tmp1 = cell2mat(cellfun(@(x) x.wt_uncon(2,:)', sim, 'UniformOutput', false));
        tmp_corr = NaN(spec.H+1,1);
        for hh = 1:spec.H+1
            tmp_corr(hh) = corr(tmp(hh,:)',tmp1(hh,:)');
        end
        plot(0:spec.H, tmp_corr, 'color', 0.6*[1,1,1], 'linewidth', lw);
        title({'Correlation Between','Cond. and Uncond. Weights'})
    elseif ii == 3
        for mm = 1:3
            if mm < 3
                col_tmp = col(mm,:); lt = '-';
            else
                col_tmp = [0,0,0]; lt = '--';
            end
            irf_mean = mean(cell2mat(cellfun(@(x) x.irf_mean{mm}, sim, 'UniformOutput', false)),2,'omitnan');
            plot(0:spec.H,irf_mean - sim_par.rho_y.^(0:20)',lt,'linewidth',lw,'color',col_tmp)
            hold on
        end
        irf_mean = mean(cell2mat(cellfun(@(x) x.irf_mean_uncon, sim, 'UniformOutput', false)),2,'omitnan');
        plot(0:spec.H,irf_mean - sim_par.rho_y.^(0:20)','--','linewidth',lw,'color',col(end-1,:))
        title('Average Bias')
        xlabel('Horizon')
    elseif ii == 4
        for mm = 1:3
            if mm < 3
                col_tmp = col(mm,:); lt = '-';
            else
                col_tmp = [0,0,0]; lt = '--';
            end
            irf_std = mean(cell2mat(cellfun(@(x) x.irf_std{mm}, sim, 'UniformOutput', false)),2,'omitnan');
            plot(0:spec.H,irf_std,lt,'linewidth',lw,'color',col_tmp)
            hold on
        end
        irf_std = mean(cell2mat(cellfun(@(x) x.irf_std_uncon, sim, 'UniformOutput', false)),2,'omitnan');
        plot(0:spec.H,irf_std,'--','linewidth',lw,'color',col(end-1,:))
        title('Average Standard Deviation')
        xlabel('Horizon')
        legend('VAR','LP','Average (Cond.)', 'Average (Uncond.)',...
            'position',[0.79,0.11,0.2,0.2],'fontsize',11-sl)
    end
    hold off
    set(gca,'fontsize',12-sl)
    xlim([0,20])
    grid on

end

% set(ff,'Units','Inches');
% pos = get(ff,'Position');
% set(ff,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% if sl == 1
%     print(ff, 'figures/sim_ARMA_slides','-dpdf','-bestfit')
% else
%     print(ff, 'figures/sim_ARMA','-dpdf','-bestfit')
% end
