% ========================================================================
% PLOTS FOR MONETARY POLICY SHOCK EMPIRICAL APPLICATION (SECTION 4.1 AND 
% APPENDIX D)
% ========================================================================

clear
clc
addpath(genpath(pwd))

%%

% ------------------------------------------------------------------------
% BASELINE (SECTION 4.1)
% ------------------------------------------------------------------------

clear
clc
addpath(genpath(pwd))
load('out/Ramey_Monthlydat6996_2023-08-01.mat')

% PLOT COLORS
col = [76, 114, 176;
    221, 132, 82;
    85, 168, 104;
    129, 114, 179;
    196, 78, 82;
    204, 185, 116;
    218, 139, 195;
    100,181,205]/255;

% PRELIMINARIES
H = spec.H;
sc = [100, 1, 100, 100]/4;
var_names = {'IP', 'Unemployment', 'CPI', 'CPI (Commodities)'};
mod_names = {'Weighted Average', 'VAR', ...
    'LP (Contemp. Controls)', 'LP (Only Lagged Controls)'};

% PLOT FIGURE 3
fs = 12;
ff = figure('position',[50,50,900,900]);
for ii = 1:3
    
    % WEIGHTS
    subplot('position', [0.09,1.01-ii*0.15,0.25,0.11])
    hold on
    for mm = 3:-1:1
        plot(0:H, squeeze(wt.wt_blk(mm,ii,:)), 'color', col(mm,:), 'linewidth', 1.25)
    end
    hold off
    xlim([0,H])
    grid on
    ylabel(var_names{ii},'FontWeight','bold','FontSize',fs+2)
    set(gca,'fontsize',fs)
    if ii == 1
        title('Weights','FontSize',fs+1)
    elseif ii == 3
        xlabel('Horizon (Months)','Fontsize',fs)
    end
    xticks(0:10:50)
    box on
end

% IMPULSE RESPONSES
for mm = 1:3
    
    for ii = 1:3

        rr = floor(mm/2);
        cc = mod(mm,2);
    
        subplot('position', [0.09+cc*0.32, 1.01-ii*0.15-rr*0.5, 0.25, 0.11])
        
        hold on

        % model IRF
        irf_tmp = sc(ii)*squeeze(est{mm}.irf{end}(ii,1,:,:));
        irf_m = mean(irf_tmp,2);
        irf_q = quantile(irf_tmp,normcdf([-1,1]),2);
    
        % average IRF
        irf_av = sc(ii)*squeeze(wt.irf_av(ii,:,:));
        irf_av_m = mean(irf_av,2);
        irf_av_q = quantile(irf_av,normcdf([-1,1]),2);
        p(1) = plot(0:H, mean(irf_av,2), '-', 'Color', 'k', 'linewidth', 1.5);
        hold on
        fill([0:H, H:-1:0],[irf_av_q(:,2)', fliplr(irf_av_q(:,1)')],...
            'k', 'FaceAlpha', 0.15, 'LineWidth', 0.25, 'EdgeAlpha', 0.5, 'LineStyle', 'none')
        p(mm+1) = plot(0:H, irf_m, '-', 'Color', col(mm,:), 'LineWidth', 1.5);
        plot(0:H, irf_q, '--', 'Color', col(mm,:), 'LineWidth', 0.5);

        if ii == 1
            ylim([-1,0.5])
        elseif ii == 2
            ylim([-0.2,0.3])
        elseif ii == 3
            ylim([-1,0.3])
        elseif ii == 4
            ylim([-10,3.5])
        end
        xlim([0,H])
        xticks(0:10:50)
        grid on
        hold off
        set(gca,'fontsize',fs)
        if ii == 1
            title(strcat('IRF:',{' '},mod_names{mm+1}),'FontSize',fs+1)
        elseif ii == 3
            xlabel('Horizon (Months)','Fontsize',fs)
        end
        box on

        if mm == 2
            ylabel(var_names{ii},'FontWeight','bold','FontSize',fs+2)
        end
        if mm == 3 && ii == 1
            legend(p,mod_names,...
                'position', [0.69, 0.06, 0.3, 0.15], 'fontsize', 11)
        end
    end
end

% set(ff,'Units','Inches');
% pos = get(ff,'Position');
% set(ff,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(ff, 'figures/Ramey_baseline','-dpdf','-bestfit')


%%

% ------------------------------------------------------------------------
% DIFFERENT SAMPLE SPLITS (APPENDIX D.1)
% ------------------------------------------------------------------------

% LOAD OUTPUT
clear
out_0 = load('out/Ramey_Monthlydat6996-0blks_2023-08-01.mat');
out_2 = load('out/Ramey_Monthlydat6996_2023-08-01.mat');
out_5 = load('out/Ramey_Monthlydat6996-5blks_2023-08-01.mat');

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

% PRELIMINARIES
H = out_2.spec.H;
sc = [100, 1, 100, 100];
var_names = {'IP', 'Unemployment', 'CPI', 'CPI (Commodities)'};

% PLOT FIGURE A.4 (PREDICTION POOL WEIGHTS)
ff = figure('position',[50,50,900,575]);
for ii = 1:3

    % SELECT INSTRUMENT
    if ii == 1
        wt = out_0.wt.wt_glb;
        blk_name = 'In-sample';
    elseif ii == 2
        wt = out_2.wt.wt_blk;
        blk_name = '2 Blocks';
    elseif ii == 3
        wt = out_5.wt.wt_blk;
        blk_name = '5 Blocks';
    end

    % PLOT
    for jj = 1:3
        subplot('position', [ii*0.23-0.16, 0.99-jj*0.3, 0.19, 0.21])
        colororder(col);
        plot(0:H, squeeze(wt(:,jj,:))','linewidth',1.25)
        xlim([0,H])
        xticks(0:10:50)
        grid on
        if ii == 1
            ylabel(var_names{jj},'fontweight','bold')
        end
        set(gca,'fontsize',11)
        if jj == 1
            title(blk_name)
            if ii == 3
                legend('VAR','LP (Contemp. Controls)','LP (Only Lagged Controls)',...
                    'position', [0.75,0.09,0.23,0.17],'fontsize',11)
            end
        end
        if jj == 3
            xlabel('Horizon (Months)')
        end
    end

end

% set(ff,'Units','Inches');
% pos = get(ff,'Position');
% set(ff,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(ff, 'figures/Ramey_split_weights','-dpdf','-bestfit')

% PLOT FIGURE A.5 (AVERAGE IMPULSE RESPONSES)
ff = figure('position',[50,50,650,550]);
for ii = 1:2

    % SELECT INSTRUMENT
    if ii == 1
        alt = out_0;
        blk_name = 'In-sample';
    elseif ii == 2
        alt = out_5;
        blk_name = '5 Blocks';
    end

    % PLOT
    for jj = 1:3
        subplot('position', [ii*0.35-0.26, 1.05-jj*0.32, 0.28, 0.23])
        colororder(col);
        hold on
        for mm = 1:2
            if mm == 1
                tmp = out_2;
                col_tmp = 0*[1,1,1];
            else
                tmp = alt;
                col_tmp = col(8+ii,:);
            end
            irf_tmp = sc(jj)*squeeze(tmp.wt.irf_av(jj,:,:));
            lw = 1.5; lt = '-';
            irf_m = mean(irf_tmp,2);
            irf_q = quantile(irf_tmp,normcdf([-1,1]),2);
    
            plt = plot(0:H, irf_m, lt, 'color', col_tmp, 'linewidth', lw);
            if ii == 1
                if mm == 2
                    p(1) = plt;
                elseif mm == 1
                    p(2) = plt;
                end
            elseif ii == 2 && mm == 2
                p(3) = plt;
            end
            fill([0:H, H:-1:0],[irf_q(:,2)', fliplr(irf_q(:,1)')],...
                col_tmp, 'FaceAlpha', 0.15,...
                'LineWidth', 0.25, 'EdgeAlpha', 0.5, 'LineStyle', 'none')
        end
        box on
        if jj == 1
            ylim([-3.2,1.5])
        elseif jj == 2
            ylim([-0.5,0.8])
        elseif jj == 3
            ylim([-2.7,1])
        elseif jj == 4
            ylim([-100,50])
        end
        xlim([0,H])
        xticks(0:10:50)
        grid on
        if ii == 1
            ylabel(var_names{jj},'fontweight','bold')
        end
        set(gca,'fontsize',11)
        if jj == 1
            title(blk_name)
            if ii == 2
                legend(p, {'In-sample', '2 Blocks', '5 Blocks'},...
                    'position', [0.75,0.09,0.23,0.15],'fontsize',11)
            end
        end
        if jj == 3
            xlabel('Horizon (Months)')
        end
    end

end

% set(ff,'Units','Inches');
% pos = get(ff,'Position');
% set(ff,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(ff, 'figures/Ramey_split_IRFs','-dpdf','-bestfit')


%%

% ------------------------------------------------------------------------
% INCLUDE CHOLESKY VAR INTO POOL (APPENDIX D.2)
% ------------------------------------------------------------------------

% LOAD OUTPUT
clear
load('out/Ramey_Monthlydat6996-4mods_2023-08-08.mat')
bm = load('out/Ramey_Monthlydat6996_2023-08-01.mat');

% PLOT COLORS
col = [218, 139, 195;
    76, 114, 176;
    221, 132, 82;
    85, 168, 104;
    129, 114, 179;
    196, 78, 82;
    204, 185, 116;
    100,181,205;
    0,120,140]/255;

% PRELIMINARIES
H = spec.H;
sc = [100, 1, 100, 100]/4;
var_names = {'IP', 'Unemployment', 'CPI', 'CPI (Commodities)'};
mod_names = {'Average (3 Models)', ...
    'Average (4 Models)', ...
    'VAR (Cholesky)', ...
    'VAR (Int Instrument)',...
    'LP (Contemp. Controls)', 'LP (Only Lagged Controls)'};

% PLOT FIGURE A.6
fs = 12;
ff = figure('position',[50,50,800,900]);
for ii = 1:3
    
    % WEIGHTS
    subplot('position', [0.1,1.01-ii*0.15,0.25,0.11])
    hold on
    for mm = 4:-1:1
        plot(0:H, squeeze(wt.wt_blk(mm,ii,:)), 'color', col(mm,:), 'linewidth', 1.25)
    end
    hold off
    xlim([0,H])
    grid on
    ylabel(var_names{ii},'FontWeight','bold','FontSize',fs+2)
    set(gca,'fontsize',fs)
    if ii == 1
        title('Weights','FontSize',fs+1)
    elseif ii == 3
        xlabel('Horizon (Months)','Fontsize',fs)
    end
    xticks(0:10:50)
    box on
end

for ii = 1:3

    % average IRF
    irf_av = sc(ii)*squeeze(wt.irf_av(ii,:,:));
    irf_av_m = mean(irf_av,2);
    irf_av_q = quantile(irf_av,normcdf([-1,1]),2);
    irf_av_bm = sc(ii)*squeeze(bm.wt.irf_av(ii,:,:));
    irf_av_bm_m = mean(irf_av_bm,2);
    irf_av_bm_q = quantile(irf_av_bm,normcdf([-1,1]),2);

    % PLOT
    cc = 1; rr = 0;
    subplot('position', [0.1+cc*0.33, 1.01-ii*0.15-rr*0.5, 0.25, 0.11])

    p(1) = plot(0:H, irf_av_bm_m, '-', 'Color', 'k', 'linewidth', 1.5);
    hold on
    p(2) = plot(0:H, irf_av_m, '-', 'Color', col(9,:), 'linewidth', 1.5);
    fill([0:H, H:-1:0],[irf_av_bm_q(:,2)', fliplr(irf_av_bm_q(:,1)')],...
        'k', 'FaceAlpha', 0.15, 'LineWidth', 0.25, 'EdgeAlpha', 0.5, 'LineStyle', 'none')
    fill([0:H, H:-1:0],[irf_av_q(:,2)', fliplr(irf_av_q(:,1)')],...
        col(9,:), 'FaceAlpha', 0.15, 'LineWidth', 0.25, 'EdgeAlpha', 0.5, 'LineStyle', 'none')
    
    % SPECS
    if ii == 1
            ylim([-1,0.5])
        elseif ii == 2
            ylim([-0.2,0.3])
        elseif ii == 3
            ylim([-1,0.3])
        elseif ii == 4
            ylim([-10,3.5])
    end
    xlim([0,H])
    xticks(0:10:50)
    grid on
    hold off
    set(gca,'fontsize',fs)
    if ii == 1
        title('IRF: Weighted Average','FontSize',fs+1)
    elseif ii == 3
        xlabel('Horizon (Months)','Fontsize',fs)
    end
    box on

end

% IMPULSE RESPONSES
for mm = 1:2
    
    for ii = 1:3

        rr = 1; cc = mod(mm+1,2);
    
        subplot('position', [0.09+cc*0.33, 1.01-ii*0.15-rr*0.5, 0.25, 0.11])
        
        hold on

        % model IRF
        irf_1 = sc(ii)*squeeze(est{2*(mm-1)+1}.irf{end}(ii,1,:,:));
        irf_m_1 = mean(irf_1,2);
        irf_q_1 = quantile(irf_1,normcdf([-1,1]),2);

        irf_2 = sc(ii)*squeeze(est{2*(mm-1)+2}.irf{end}(ii,1,:,:));
        irf_m_2 = mean(irf_2,2);
        irf_q_2 = quantile(irf_2,normcdf([-1,1]),2);

        p(mm*2+1) = plot(0:H, irf_m_1, '-', 'Color', col(2*(mm-1)+1,:), 'LineWidth', 1.5);
        p(mm*2+2) = plot(0:H, irf_m_2, '-', 'Color', col(2*(mm-1)+2,:), 'LineWidth', 1.5);
        plot(0:H, irf_q_1, '--', 'Color', col(2*(mm-1)+1,:), 'LineWidth', 0.5);
        plot(0:H, irf_q_2, '--', 'Color', col(2*(mm-1)+2,:), 'LineWidth', 0.5);

        if ii == 1
            ylim([-1,0.5])
        elseif ii == 2
            ylim([-0.2,0.3])
        elseif ii == 3
            ylim([-1,0.3])
        elseif ii == 4
            ylim([-10,3.5])
        end
        xlim([0,H])
        xticks(0:10:50)
        grid on
        hold off
        set(gca,'fontsize',fs)
        if ii == 1
            if mm == 1
                title('IRF: VARs','FontSize',fs+1)
            elseif mm == 2
                title('IRF: LPs','FontSize',fs+1)
            end
        elseif ii == 3
            xlabel('Horizon (Months)','Fontsize',fs)
        end
        box on

        if mm == 1
            ylabel(var_names{ii},'FontWeight','bold','FontSize',fs+2)
        end
        if mm == 2 && ii == 1
            legend(p,mod_names,...
                'position', [0.69, 0.06, 0.3, 0.25], 'fontsize', 11)
        end
    end
end

% set(ff,'Units','Inches');
% pos = get(ff,'Position');
% set(ff,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(ff, 'figures/Ramey_4mods','-dpdf','-bestfit')

%%

% ------------------------------------------------------------------------
% ALTERNATIVE INSTRUMENTS (APPENDIX D.3)
% ------------------------------------------------------------------------

% LOAD OUTPUT
clear
JK  = load('out/Ramey_JK_2023-08-08.mat');
GK  = load('out/Ramey_GK_2023-08-01.mat');
MAR = load('out/Ramey_MAR_2023-08-08.mat');

% PLOT COLORS
col = [76, 114, 176;
    221, 132, 82;
    85, 168, 104;
    129, 114, 179;
    196, 78, 82;
    204, 185, 116;
    218, 139, 195;
    100,181,205]/255;

% PRELIMINARIES
H = JK.spec.H;
sc = [100, 1, 100, 100]/4;
var_names = {'IP', 'Unemployment', 'CPI', 'CPI (Commodities)'};


% PLOT FIGURE A.7 (OPTIMAL WEIGHTS)
ff = figure('position',[50,50,900,575]);
for ii = 1:3

    % SELECT INSTRUMENT
    if ii == 1
        tmp = GK;
        inst_name = 'Gertler and Karadi';
    elseif ii == 2
        tmp = MAR;
        inst_name = {'Miranda-Agrippino','and Ricco'};
    elseif ii == 3
        tmp = JK;
        inst_name = 'Jarociński and Karadi';
    end

    % PLOT
    for jj = 1:3
        subplot('position', [ii*0.23-0.16, 0.99-jj*0.3, 0.19, 0.21])
        colororder(col);
        plot(0:H, squeeze(tmp.wt.wt_blk(:,jj,:))','linewidth',1.25)
        xlim([0,H])
        xticks(0:10:50)
        grid on
        if ii == 1
            ylabel(var_names{jj},'fontweight','bold')
        end
        if jj == 1
            title(inst_name)
            legend('VAR','LP (Contemp. Controls)','LP (Only Lagged Controls)',...
                'position', [0.75,0.09,0.23,0.15],'fontsize',11)
        elseif jj == 3
            xlabel('Horizon (Months)')
        end
        set(gca,'fontsize',11)
    end

end

% set(ff,'Units','Inches');
% pos = get(ff,'Position');
% set(ff,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(ff, 'figures/Ramey_alt_weights','-dpdf','-bestfit')

% PLOT FIGURE A.8 (AVERAGE IMPULSE RESPONSES)
ff = figure('position',[50,50,900,575]);
for ii = 1:3

    % SELECT INSTRUMENT
    if ii == 1
        tmp = GK;
        inst_name = 'Gertler and Karadi';
    elseif ii == 2
        tmp = MAR;
        inst_name = {'Miranda-Agrippino','and Ricco'};
    elseif ii == 3
        tmp = JK;
        inst_name = 'Jarociński and Karadi';
    end

    % PLOT
    for jj = 1:3
        subplot('position', [ii*0.23-0.16, 0.99-jj*0.3, 0.19, 0.21])
        colororder(col);
        hold on
        for mm = 3:-1:0
            if mm > 0
                irf_tmp = sc(jj)*squeeze(tmp.est{mm}.irf{end}(jj,1,:,:));
                col_tmp = col(mm,:); lw = 1.25; lt = '-';
            else
                irf_tmp = sc(jj)*squeeze(tmp.wt.irf_av(jj,:,:));
                col_tmp = 0*[1,1,1]; lw = 1.5; lt = '--';
            end
            irf_m = mean(irf_tmp,2);
            irf_q = quantile(irf_tmp,normcdf([-1,1]),2);
    
            p(mm+1) = plot(0:H, irf_m, lt, 'color', col_tmp, 'linewidth', lw);
            if mm > 0
                plot(0:H, irf_q, ':', 'Color', col_tmp, 'LineWidth', 0.5)
            else
                fill([0:H, H:-1:0],[irf_q(:,2)', fliplr(irf_q(:,1)')],...
                    col_tmp, 'FaceAlpha', 0.15,...
                    'LineWidth', 0.25, 'EdgeAlpha', 0.5, 'LineStyle', 'none')
            end
        end
        if jj == 1
            ylim([-2.5,1.5])
        elseif jj == 2
            ylim([-1,1])
            yticks(-1:0.5:1)
        elseif jj == 3
            ylim([-1.0,0.2])
        end
        xlim([0,H])
        xticks(0:10:50)
        grid on
        box on
        if ii == 1
            ylabel(var_names{jj},'fontweight','bold')
        end
        if jj == 1
            title(inst_name)
            legend(p, {'Weighted Average','VAR','LP (Contemp. Controls)','LP (Only Lagged Controls)'},...
                'position', [0.75,0.09,0.23,0.2],'fontsize',11)
        elseif jj == 3
            xlabel('Horizon (Months)')
        end
        set(gca,'fontsize',11)
    end

end

% set(ff,'Units','Inches');
% pos = get(ff,'Position');
% set(ff,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(ff, 'figures/Ramey_alt_IRFs','-dpdf','-bestfit')