% ========================================================================
% PLOTS FOR MONTE CARLO WITH AR(1) MODEL (SECTION 3.1 AND APPENDIX B)
% ========================================================================

clear
clc
addpath(genpath(pwd))

%%

% ------------------------------------------------------------------------
% BASELINE (SECTION 3.1)
% ------------------------------------------------------------------------

% LOAD OUTPUT
clear
load('out/sim_HJ_lean_2023-09-21.mat')

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

% PLOT FIGURE 1
ff = figure('position',[100,100,700,500]);
lw = 1.5;
for ii = 1:4
    colororder(col)
    if ii ~= 2
        subplot('position', [mod(ii+1,2)*0.49+0.07, 1.07-ceil(ii/2)*0.48, 0.4, 0.35])
    end
    if ii == 1
        wt    = mean(cell2mat(cellfun(@(x) x.wt(1,:), sim, 'UniformOutput', false)));
        wt_LS = mean(cell2mat(cellfun(@(x) x.wt_LS(1,:), sim, 'UniformOutput', false)));
        plot(0:spec.H,wt,'k--','linewidth',lw)
        hold on
        plot(0:spec.H,wt_LS,'--','color',col(end-1,:),'linewidth',lw)
        title('Weight on Model With Controls')
        ylim([0,1])
    elseif ii == 3
        for mm = 1:spec.n_model+2
            if mm <= spec.n_model
                col_tmp = col(mm,:); lt_tmp = '-';
            elseif mm == 3; col_tmp = [0,0,0]; lt_tmp = '--';
            elseif mm == 4; col_tmp = col(end-1,:); lt_tmp = '--';
            end
            irf_mean = mean(cell2mat(cellfun(@(x) x.irf_mean{mm}', sim, 'UniformOutput', false)),'omitnan');
            plot(0:spec.H,irf_mean - sim_par.rho.^(0:20),lt_tmp,'linewidth',lw,'color',col_tmp)
            hold on
        end
        title('Average Bias')
        xlabel('Horizon')
    elseif ii == 4
        for mm = 1:spec.n_model+2
            if mm <= spec.n_model
                col_tmp = col(mm,:); lt_tmp = '-';
            elseif mm == 3; col_tmp = [0,0,0]; lt_tmp = '--';
            elseif mm == 4; col_tmp = col(end-1,:); lt_tmp = '--';
            end
            irf_std = std(cell2mat(cellfun(@(x) x.irf_mean{mm}', sim, 'UniformOutput', false)),'omitnan');
            plot(0:spec.H,irf_std,lt_tmp,'linewidth',lw,'color',col_tmp)
            hold on
        end
        title('Std. Dev. of Point Estimates')
        xlabel('Horizon')
        legend('With Controls','No Controls', 'Prediction Pool Average', 'Least Squares Average',...
            'position',[0.56,0.59,0.32,0.2],'fontsize',11)
    end
    hold off
    set(gca,'fontsize',12)
    xlim([0,20])
    grid on

end

% set(ff,'Units','Inches');
% pos = get(ff,'Position');
% set(ff,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(ff, 'figures/sim_HJ','-dpdf','-bestfit')


%%

% ------------------------------------------------------------------------
% POOL WITH VAR (APPENDIX B)
% ------------------------------------------------------------------------

% LOAD OUTPUT
clear
load('out/sim_HJ_VAR_lean_2023-09-25.mat')

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

% PLOT FIGURE A.1
ff = figure('position',[100,100,700,500]);
lw = 1.5;
for ii = 1:4
    colororder(col)
    if ii~=2
        subplot('position', [mod(ii+1,2)*0.49+0.07, 1.07-ceil(ii/2)*0.48, 0.4, 0.35])
    end
    if ii == 1
        wt = reshape(cell2mat(cellfun(@(x) x.wt, sim, 'UniformOutput', false)),3,length(sim),spec.H+1);
        wt = squeeze(mean(wt,2));
        plot(0:spec.H,wt,'linewidth',lw)
        title('Weights')
        ylim([0,1])
    elseif ii == 3
        for mm = 1:spec.n_model+1
            if mm <= spec.n_model
                col_tmp = col(mm,:); lt_tmp = '-';
            elseif mm == 4; col_tmp = [0,0,0]; lt_tmp = '--';
            end
            irf_mean = mean(cell2mat(cellfun(@(x) x.irf_mean{mm}', sim, 'UniformOutput', false)),'omitnan');
            plot(0:spec.H,irf_mean - sim_par.rho.^(0:20),lt_tmp,'linewidth',lw,'color',col_tmp)
            hold on
        end
        title('Average Bias')
    elseif ii == 4
        for mm = 1:spec.n_model+1
            if mm <= spec.n_model
                col_tmp = col(mm,:); lt_tmp = '-';
            elseif mm == 4; col_tmp = [0,0,0]; lt_tmp = '--';
            end
            irf_std = std(cell2mat(cellfun(@(x) x.irf_mean{mm}', sim, 'UniformOutput', false)),'omitnan');
            plot(0:spec.H,irf_std,lt_tmp,'linewidth',lw,'color',col_tmp)
            hold on
        end
        legend('LP: With Controls', 'LP: No Controls', 'VAR', 'Average',...
            'position',[0.56,0.59,0.27,0.2],'fontsize',11)
        title('Std. Dev. of Point Estimates')
        xlabel('Horizon')
    end
    hold off
    set(gca,'fontsize',12)
    xlim([0,20])
    grid on

end

% set(ff,'Units','Inches');
% pos = get(ff,'Position');
% set(ff,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(ff, 'figures/sim_HJ_VAR','-dpdf','-bestfit')

