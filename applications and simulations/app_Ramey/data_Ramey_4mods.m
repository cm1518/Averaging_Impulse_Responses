function data = data_Ramey_4mods(file_name, sheet_name, shock)

% ========================================================================
%
% READ IN DATA: MONETARY SHOCK APPLICATION WITH CHOLESKY VAR
%
% ========================================================================

% READ DATA
data.tab = readtable(file_name, 'Sheet', sheet_name);
    
% FULL LIST OF VARIABLE ABBREVIATIONS
data.var_list_all = {'LIP', 'UNEMP', 'LCPI', 'LPCOM', 'FFR', strcat('CUM',shock) ,shock};

% FULL LIST OF LONG VARIABLES NAMES
data.var_list_long_all = [...
    "Log Industrial Production", ...
    "Unemployment", ...
    "Log Consumer Price Index", ...
    "Log Commodity Price Index",...
    "Federal Funds Rate",...
    "Cumulation of Shock",...
    "Shock"];

% REMOVE NaN ROWS
data.tab = data.tab(:,ismember(data.tab.Properties.VariableNames,data.var_list_all));
data.tab = data.tab(~any(ismissing(data.tab),2),:);

% NUMBER OF VARIABLES
data.n_var_total = size(data.var_list_all,2);

% NUMBER OF PERIODS
data.T = size(data.tab,1);

end