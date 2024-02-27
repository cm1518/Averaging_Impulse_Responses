% ========================================================================
% PLOTS FOR TFP SHOCK EMPIRICAL APPLICATION (SECTION 4.2)
% ========================================================================

clear
clc
addpath(genpath(pwd))

% LOAD OUTPUT
load('out/TFP_2023-08-03.mat')

% PLOT COLORS
col = [76, 114, 176;
    221, 132, 82;
    85, 168, 104;
    129, 114, 179;
    196, 78, 82;
    204, 185, 116;
    218, 139, 195;
    100,181,205;
    0,120,140]/255;

% PRELIMINARIES
H = spec.H;
model_names = {'VAR (All Controls)',...
    'VAR (Factors)',...
    'VAR (No Controls)',...
    'FAVAR',...
    'LP (All Controls)',...
    'LP (Factors)',...
    'LP (No Controls)',...
    'Single Equation'};

% PLOT FIGURE 4 (IRFs)
ff = figure('position',[100,100,700,700]);
for mm = 1:length(est)

    subplot('position', [mod(mm+1,2)*0.35+0.05, 1.04-ceil(mm/2)*0.24, 0.28, 0.16])
    
    % average
    p(1) = plot(0:H, quantile(squeeze(wt.irf_av),normcdf(0),2),'-','color',0*[1,1,1],'linewidth',1.5);
    hold on
    irf_q = quantile(squeeze(wt.irf_av),normcdf([-1,1]),2);
    fill([0:H, H:-1:0],[irf_q(:,2)', fliplr(irf_q(:,1)')],...
        0*[1,1,1], 'FaceAlpha', 0.15,...
        'LineWidth', 0.25, 'EdgeAlpha', 0.5, 'LineStyle', 'none')
    
    % model
    p(mm+1) = plot(0:H, quantile(squeeze(est{mm}.irf{3}(1,1,:,:)),normcdf(0),2),'-','color',col(mm,:),'linewidth',1.5);
    plot(0:H, quantile(squeeze(est{mm}.irf{3}(1,1,:,:)),normcdf([-1,1]),2),'--','color',col(mm,:),'linewidth',0.75)
    
    xlim([0,H])
    ylim([-0.2,0.12])
    grid on
    title(model_names{mm})
    set(gca,'fontsize',11)
    xticks(0:5:20)

    if mm == 7 || mm == 8
        xlabel('Horizon (Quarters)')
    end
    if mm == 8
        legend(p, {'Average',model_names{1:8}},'position',[0.7,0.08,0.295,0.32],'fontsize',10)
    end

end

% set(ff,'Units','Inches');
% pos = get(ff,'Position');
% set(ff,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(ff, 'figures/TFP_IRFs','-dpdf','-bestfit')

% PLOT FIGURE 5 (CONDITIONAL VS UNCONDITIONAL WEIGHTS AND IRFs)
ff = figure('position',[100,100,600,800]);
for ii = 1:3
    subplot('position',[0.08,1.035-ii*0.32,0.5,0.24])
    if ii == 1
        
        % conditional
        p_con(1) = plot(0:H, quantile(squeeze(wt.irf_av),normcdf(0),2),'-','color',0*[1,1,1],'linewidth',1.5);
        hold on
        irf_q = quantile(squeeze(wt.irf_av),normcdf([-1,1]),2);
        fill([0:H, H:-1:0],[irf_q(:,2)', fliplr(irf_q(:,1)')],...
            0*[1,1,1], 'FaceAlpha', 0.1,...
            'LineWidth', 0.25, 'EdgeAlpha', 0.5, 'LineStyle', 'none')

        % unconditional
        p_con(2) = plot(0:H, quantile(squeeze(wt.irf_av_uncon),normcdf(0),2),'-','color',col(9,:),'linewidth',1.5);
        hold on
        irf_q = quantile(squeeze(wt.irf_av_uncon),normcdf([-1,1]),2);
        fill([0:H, H:-1:0],[irf_q(:,2)', fliplr(irf_q(:,1)')],...
            col(9,:), 'FaceAlpha', 0.1,...
            'LineWidth', 0.25, 'EdgeAlpha', 0.5, 'LineStyle', 'none')
        ylim([-0.07,0.1])
        title('IRF: Conditional vs Unconditional')
        legend(p_con,{'Conditional on Shock','Unconditional'},'position',[0.61,0.715,0.36,0.075],'fontsize',10)
        
    elseif ii == 2
        plot(0:H, squeeze(wt.wt_blk),'linewidth',1.25)
        colororder(col)
        title('Weights: Conditional on Shock')
    elseif ii == 3
        plot(0:H, squeeze(wt.wt_uncon_blk),'linewidth',1.25)
        title('Weights: Unconditional')
        legend(model_names,'position',[0.61,0.075,0.37,0.26],'fontsize',10)
        xlabel('Horizon (Quarters)')
    end
    grid on
    set(gca,'fontsize',11)
end

% set(ff,'Units','Inches');
% pos = get(ff,'Position');
% set(ff,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(ff, 'figures/TFP_weights','-dpdf','-bestfit')
