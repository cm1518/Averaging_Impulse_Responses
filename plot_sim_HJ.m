% ========================================================================
% PLOTS FOR MONTE CARLO WITH AR(1) MODEL (SECTION 3.1 AND APPENDIX B)
% ========================================================================

clear
clc
addpath(genpath(pwd))

%%

% ------------------------------------------------------------------------
% BASELINE (SECTION 3.1) AND OMITTED VARIABLES (APPENDIX B.2)
% ------------------------------------------------------------------------

% LOAD OUTPUT
clear
load('out/sim_HJ_lean_2024-02-10.mat')  % Figure 1
% load('out/sim_HJ_OV_neg_lean_2024-02-04.mat') % Figure A.2
% load('out/sim_HJ_OV_pos_lean_2024-02-06.mat') % Figure A.3

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

% EXTRACT QUANTILES FROM CELLS
qtl_1 = cell2mat(cellfun(@(x) x.irf_qtl{1}', sim, 'UniformOutput', false));
qtl_2 = cell2mat(cellfun(@(x) x.irf_qtl{2}', sim, 'UniformOutput', false));
qtl_3 = cell2mat(cellfun(@(x) x.irf_qtl{3}', sim, 'UniformOutput', false));
qtl_4 = cell2mat(cellfun(@(x) x.irf_qtl{4}', sim, 'UniformOutput', false));

% PLOT FIGURE 1/A.2/A.3
eb = [0.68, 0.95];
ff = figure('position',[100,50,700,750]);
lw = 1.5;
for ii = 1:6
    colororder(col)
    if ii ~= 2
        subplot('position', [mod(ii+1,2)*0.49+0.07, 1.04-ceil(ii/2)*0.32, 0.4, 0.24])
    end
    if ii == 1
        % WEIGHTS
        wt    = mean(cell2mat(cellfun(@(x) x.wt(1,:), sim, 'UniformOutput', false)));
        wt_LS = mean(cell2mat(cellfun(@(x) x.wt_LS(1,:), sim, 'UniformOutput', false)));
        plot(0:spec.H,wt,'k--','linewidth',lw)
        hold on
        plot(0:spec.H,wt_LS,'--','color',col(end-1,:),'linewidth',lw)
        title('Weight on Model With Controls')
        ylim([0,1])
    elseif ii == 3
        % BIAS
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
    elseif ii == 4
        % STANDARD DEVIATION
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
        legend('With Controls','No Controls', 'Prediction Pool Average', 'Least Squares Average',...
            'position',[0.56,0.72,0.32,0.16],'fontsize',11)
    elseif ii == 5 || ii == 6
        % COVERAGE
        qq = norminv(1-(1-eb(ii-4))/2);
        plot(0:spec.H,mean(qtl_1<normcdf(qq) & qtl_1>normcdf(-qq)),'color',col(1,:),'linewidth',lw)
        hold on
        plot(0:spec.H,mean(qtl_2<normcdf(qq) & qtl_2>normcdf(-qq)),'color',col(2,:),'linewidth',lw)
        plot(0:spec.H,mean(qtl_3<normcdf(qq) & qtl_3>normcdf(-qq)),'k--','linewidth',lw)
        plot(0:spec.H,mean(qtl_4<normcdf(qq) & qtl_4>normcdf(-qq)),'--','color',col(end-1,:),'linewidth',lw)
        yline(normcdf(qq)-normcdf(-qq),'-','color',0.4*[1,1,1],'linewidth',1)
        hold off
        grid on
        xlabel('Horizon')
        if ii == 5
            title('Coverage (68%)')
        elseif ii == 6
            title('Coverage (95%)')
        end
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
% POOL WITH VAR (APPENDIX B.1)
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

