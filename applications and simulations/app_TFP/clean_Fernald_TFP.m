function [dTFP, cal] = clean_Fernald_TFP(file_name, smpl)

TFP = readtable(file_name,'Sheet','quarterly');
row_TFP = 2:size(TFP,1)-6;

% DATES
cal_TFP = datetime(datenum(table2array(TFP(row_TFP,1)),'yyyy:QQ'),'ConvertFrom','datenum');
smpl_TFP = cal_TFP >= smpl.st & cal_TFP <= smpl.en;

% TFP DATA
dTFP = TFP.('dtfp')(row_TFP);

% SUBSAMPLE
dTFP =    dTFP(smpl_TFP);
cal  = cal_TFP(smpl_TFP);

end