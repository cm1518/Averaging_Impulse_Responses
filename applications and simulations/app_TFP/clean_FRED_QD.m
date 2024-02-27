function [X, fac, cal] = clean_FRED_QD(file_name, var_names_X, transform_X, smpl)

% READ DATA
opts = detectImportOptions(file_name);
FRED_QD = readtable(file_name,opts);
row_FRED = 3:size(FRED_QD,1); % rows with data

% TRANSFORMATIONS
dt_raw = table2array(FRED_QD(row_FRED,2:end));
tr = table2array(FRED_QD(2,2:end));
dt_tran = NaN(size(dt_raw));
dt_tran( :   ,tr==1) =               dt_raw(:,tr==1)  ;
dt_tran(2:end,tr==2) =          diff(dt_raw(:,tr==2)) ;
dt_tran(3:end,tr==3) =     diff(diff(dt_raw(:,tr==3)));
dt_tran( :   ,tr==4) =           log(dt_raw(:,tr==4)) ;
dt_tran(2:end,tr==5) =      diff(log(dt_raw(:,tr==5)));
dt_tran(3:end,tr==6) = diff(diff(log(dt_raw(:,tr==6))));
dt_tran(3:end,tr==7) = diff(dt_raw(2:end,tr==7)./dt_raw(1:end-1,tr==7) - 1);

% CALANDER AND SAMPLE
cal = table2array(FRED_QD(row_FRED,1));
smpl_FRED = cal >= smpl.st & cal <= smpl.en;
cal = cal(smpl_FRED,:);
dt_tran = dt_tran(smpl_FRED,:);

% ESTIMATE FACTORS
dt_std = dt_tran - mean(dt_tran,'omitnan');
dt_std = dt_std./std(dt_std,'omitnan');
bal = min(~isnan(dt_std));
dt_std_bal = dt_std(:,bal);
[~, fac, ~] = pca(dt_std_bal);

% EXTRACT RELEVANT VARIABLES
X = NaN(length(row_FRED),length(var_names_X));
for ii = 1:length(var_names_X)
    tmp = FRED_QD.(var_names_X{ii})(row_FRED);
    if strcmp(transform_X{ii},'growth')
        X(:,ii) = [NaN; (tmp(2:end)./tmp(1:end-1) - 1)];
    elseif strcmp(transform_X{ii},'log')
        X(:,ii) = log(tmp);
    else
        X(:,ii) = tmp;
    end
end
X = X(smpl_FRED,:);

end