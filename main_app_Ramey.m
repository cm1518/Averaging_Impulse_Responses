% ========================================================================
% EMPIRICAL APPLICATION
% SECTION 4.1: MONETARY SHOCK
% ------------------------------------------------------------------------
% - Appendix D.1: change spec.n_blk
% - Appendix D.2: add _4mods suffix to data_Ramey and spec_Ramey
% - Appendix D.3: replace shock with JK, MAR, or GK
% ========================================================================

clear
clc
addpath(genpath(pwd))

% READ DATA
shock = 'RRORIG'; %'JK'; %'MAR'; %'GK';
if strcmp(shock,'RRORIG')
    file_name = '/data/Ramey/Monetarydat.xlsx';
    sheet_name = 'Monthlydat6996';
else
    file_name = '/data/Ramey/alt_instruments.xlsx';
    sheet_name = 'data';
end
data = data_Ramey(file_name,sheet_name,shock); % data_Ramey_4mods(file_name,sheet_name,shock);

% MODEL AND ESTIMATION SPECIFICATIONS
spec_Ramey; % spec_Ramey_4mods

% VARIABLES TO COMPUTE WEIGHTS FOR
var_list = {'LIP', 'UNEMP', 'LCPI', 'LPCOM'};

% FILE NAME TO SAVE OUTPUT
if contains(file_name, 'Monetarydat')
    out_name = strcat('Ramey_',sheet_name);
elseif contains(file_name, 'alt_instruments')
    out_name = strcat('Ramey_',shock);
end
if spec.n_blk ~= 2
    out_name = strcat(out_name,'-',num2str(spec.n_blk),'blks');
end
if spec.n_model == 4
    out_name = strcat(out_name,'-4mods');
end

% MAIN COMPUTATIONS
calc_application;
