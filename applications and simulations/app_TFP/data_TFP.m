function data = data_TFP

% ========================================================================
%
% READ IN DATA: TFP APPLICATION
%
% ========================================================================

% ========================================================================
% READ AND CLEAN DATA
% ========================================================================

% SAMPLE PERIODS
smpl.st = datetime(1968, 1, 1);
smpl.en = datetime(2019,12,31);

% ------------------------------------------------------------------------
% FRED-QD DATA
% >> GDPC1: Real GDP
% >> PCECTPI: PCE Price Index
% >> TCU: Capacity Utilization: Total Industry
% >> AWHNONAG: Average Weekly Hours of Production and Nonsupervisory Emploees: Total private
% >> UNRATE: Unemployment
% ------------------------------------------------------------------------

var_names_X = {'GDPC1' ,'PCECTPI','TCU','AWHNONAG','UNRATE'};
transform_X = {'growth', 'growth',   '',     'log',      ''};
[X, fac, ~] = clean_FRED_QD('/data/TFP/FRED_QD.csv', var_names_X, transform_X, smpl);
X = X.*[100, 100, 1, 1, 1];

% ------------------------------------------------------------------------
% FERNALD TFP DATA
% ------------------------------------------------------------------------

[dTFP, ~] = clean_Fernald_TFP('/data/TFP/quarterly_tfp.xlsx', smpl);


% ========================================================================
% EXPORT DATA FOR ESTIMATION
% ========================================================================

% FULL LIST OF VARIABLE ABBREVIATIONS
data.var_list_all = {'GDP growth', 'inflation', 'utilization', 'hours', 'unemp',...
    'fac_1', 'fac_2', 'fac_3', 'TFP growth'};

% DATA
data.tab = array2table(...
    [X, fac(:,1:3), dTFP], ...
    'VariableNames', data.var_list_all);

% NUMBER OF VARIABLES
data.n_var_total = size(data.var_list_all,2);

% NUMBER OF PERIODS
data.T = size(data.tab,1);

end