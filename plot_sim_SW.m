% ========================================================================
% PLOTS FOR MONTE CARLO WITH SMETS AND WOUTERS (2007) NK MODEL
% ========================================================================

clear
clc
addpath(genpath(pwd))

% LOAD OUTPUT
load('out/sim_SW_lean_2023-10-15.mat')

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

% COMPUTE TRUE IRF
irf_Y = irf_true_SW(5,20);
irf_true = irf_Y/irf_Y(3,1);

% PLOT FIGURE A.2
lw = 1.5;
var_name = {'GDP', 'Inflation', 'Interest Rate', 'Consumption', 'Investment','Wages', 'Hours'};
ff = figure('position',[50,50,900,950]);
ll = -0.20;
tt = 1;
d_hor = 0.27;
d_ver = 0.135;
ht = 0.10;
wd = 0.21;
for ii = 1:7

    for jj = 1:3
        subplot('position', [ll+jj*d_hor, tt-d_ver*ii, wd, ht])
        if jj == 1
            for mm = 1:2
                wt = squeeze(mean(cell2mat(cellfun(@(x) x.wt(2,ii,:), sim, 'UniformOutput', false))));
                plot(1:spec.H, wt(2:end),'k--','linewidth',lw)
                hold on
            end
            hold off
            ylabel(var_name{ii},'fontweight','bold')
            ax = gca;
            ylim([0,ax.YLim(2)])
            if ii == 1
                title('Weight on LP')
            end
        elseif jj == 2
            for mm = 1:3
                if mm <= 2; col_tmp = col(mm,:); lt = '-';
                elseif mm == 3; col_tmp = [0,0,0]; lt = '--';
                end
                irf_mean = mean(cell2mat(cellfun(@(x) x.irf_mean{mm}(ii,:), sim, 'UniformOutput', false)));
                irf_bias = irf_mean - irf_true(ii,:);
                plot(1:spec.H, irf_bias(2:end), lt,'color', col_tmp,'linewidth',lw)
                hold on
            end
            hold off
            if ii == 1
                title('Average Bias')
            end
        elseif jj == 3
            for mm = 1:3
                if mm <= 2; col_tmp = col(mm,:); lt = '-';
                elseif mm == 3; col_tmp = [0,0,0]; lt = '--';
                end
                irf_std = mean(cell2mat(cellfun(@(x) x.irf_std{mm}(ii,:), sim, 'UniformOutput', false)));
                plot(1:spec.H, irf_std(2:end), lt,'color', col_tmp,'linewidth',lw)
                hold on
            end
            hold off            
            if ii == 1
                title('Average Standard Deviation')
                legend({'VAR','LP','Average'},...
                    'position', [0.84,tt-d_ver*7,0.15,0.12],'fontsize',11)
            end
        end
        
        set(gca,'fontsize',11)
        xlim([1,16])
        xticks([1,5:5:15])
        if ii == 7
            xlabel('Horizon')
        end
        grid on
    
    end

end

% set(ff,'Units','Inches');
% pos = get(ff,'Position');
% set(ff,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(ff, 'figures/sim_SW','-dpdf','-bestfit')

% PLOT FIGURE A.3
ff = figure('position',[50,50,900,350]);
for ii = 1:7
    if ii <= 3; jj = ii;
    else; jj = ii + 1;
    end
    subplot('position',[0.04+mod(jj-1,4)*0.25, 1.08-ceil(jj/4)*0.47, 0.19, 0.32])
    for mm = 1:3
        if mm <= 2; col_tmp = col(mm,:); lt = '-';
        elseif mm == 3; col_tmp = [0,0,0]; lt = '--';
        end
        irf_mean = mean(cell2mat(cellfun(@(x) x.irf_mean{mm}(ii,:), sim, 'UniformOutput', false)));
        plot(0:spec.H, irf_mean, lt, 'color' ,col_tmp, 'linewidth', lw)
        hold on
    end
    plot(0:spec.H, irf_true(ii,:),'color', 0.7*[1,1,1], 'linewidth', lw)
    hold off
    xlim([0,16])
    title(var_name{ii}) %,'fontweight','bold')
    if ii > 3        
        xlabel('Horizon')
    end
    if ii == 7
        legend({'VAR','LP','Average','True'},...
            'position', [0.04+3*0.25,1.08-0.47,0.15,0.30],'fontsize',11)
    end
    grid on
    set(gca,'fontsize',11)
end

% set(ff,'Units','Inches');
% pos = get(ff,'Position');
% set(ff,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(ff, 'figures/sim_SW_IRF','-dpdf','-bestfit')
