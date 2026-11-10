% Used for: Supplementary Tables 2-5 (grid-by-grid metrics per season and year, written to model_error_tables.xlsx) and the inputs of Extended Data Fig. 6
% Original: Code\scatterMLD_ARGOvsModels_point_season.m (data paths updated to I:\Prj; computation unchanged)
% Validate for each season per year
% Haipeng Zhao
% 4/3/2025

clear
clc

PATH_ARGO = 'I:\Prj\MLD_ARGO\one';
PATH_ECCO = 'I:\Prj\MLD_ECCO\one';
PATH_NEMO = 'I:\Prj\MLD_1M_NEMO\one_d03\mld';
PATH_HYCOM = 'I:\Prj\MLD_HYCOM\one_gat\mld';

mld_argo_array = [];
mld_nemo_array = [];
mld_ecco_array = [];
mld_hycom_array = [];

years = (2003 : 2022)'; 

%% Spring
for i = 2003 : 2022
    mld_argo = nan( 360, 180, 12, 'single' );
    for j = 3 : 5
            sub_folder1 = [PATH_ARGO, '\', '0.03sigmaTheta_ml_', num2str( i ), '.nc'];
            files = dir( [sub_folder1, '\', '*.nc'] );
            mld_n = ncread( [files(j).folder, '\', files(j).name], 'MLD' );
            mld_n(:, 1 : 90) = nan;
            mld_argo(:, :, j)  = mld_n;

            mld_s = ncread( [files(j + 6).folder, '\', files(j + 6).name], 'MLD' );
            mld_s(:, 91 : 180) = nan;
            mld_argo(:, :, j + 6)  = mld_s;
    end
    mld_argo = permute( mld_argo, [2, 1, 3] );

    d_ecco = load(  [PATH_ECCO, '\', 'mld_', num2str( i ), '.mat'] );
    d_nemo = load(  [PATH_NEMO, '\', 'mld_', num2str( i ), '.mat'] );
    d_hycom = load(  [PATH_HYCOM, '\', 'mld_', num2str( i ), '.mat'] );

    % mld_argo_array = [mld_argo_array; mld_argo(:)];
    % mld_nemo_array = [mld_nemo_array; d_nemo.mldG1deg(:)];
    % mld_ecco_array = [mld_ecco_array; d_ecco.mldG1deg(:)];
    % mld_hycom_array = [mld_hycom_array; d_hycom.mldG1deg(:)];

    mld_argo_array = [mld_argo(:)];
    mld_nemo_array = [d_nemo.mldG1deg(:)];
    mld_ecco_array = [d_ecco.mldG1deg(:)];
    mld_hycom_array = [d_hycom.mldG1deg(:)];

    % Argo vs ECCO
    idx  = ~isnan( mld_argo_array + mld_ecco_array );
    y = mld_argo_array(idx);
    x = mld_ecco_array(idx);

    x = log10( x ); 
    y = log10( y );
    N(i - 2002, 1) = length( y ); 
    mdl = fitlm( x, y );
    error_slope(i - 2002, 1) = mdl.Coefficients.Estimate(2);
    error_r2(i - 2002, 1) = mdl.Rsquared.Ordinary;
    error_bias(i - 2002, 1) = 10 .^ median( x - y ); % transformed back 
    error_mae(i - 2002, 1) = 10 .^ median( abs( x - y ) );
    error_rmse(i - 2002, 1) = 10 .^ sqrt(mean( (x - y) .^ 2 ));
    c = corrcoef( x , y );
    error_r(i - 2002, 1) = c(2);

    % Argo vs NEMO
    idx  = ~isnan( mld_argo_array + mld_nemo_array );
    y = mld_argo_array(idx);
    x = mld_nemo_array(idx);

    x = log10( x ); 
    y = log10( y );
    N(i - 2002, 2) = length( y ); 
    mdl = fitlm( x, y );
    error_slope(i - 2002, 2) = mdl.Coefficients.Estimate(2);
    error_r2(i - 2002, 2) = mdl.Rsquared.Ordinary;
    error_bias(i - 2002, 2) = 10 .^ median( x - y ); % transformed back 
    error_mae(i - 2002, 2) = 10 .^ median( abs( x - y ) );
    error_rmse(i - 2002, 2) = 10 .^ sqrt(mean( (x - y) .^ 2 ));
    c = corrcoef( x , y );
    error_r(i - 2002, 2) = c(2);

    % Argo vs HyCOM
    idx  = ~isnan( mld_argo_array + mld_hycom_array );
    y = mld_argo_array(idx);
    x = mld_hycom_array(idx);

    x = log10( x ); 
    y = log10( y );
    N(i - 2002, 3) = length( y ); 
    mdl = fitlm( x, y );
    error_slope(i - 2002, 3) = mdl.Coefficients.Estimate(2);
    error_r2(i - 2002, 3) = mdl.Rsquared.Ordinary;
    error_bias(i - 2002, 3) = 10 .^ median( x - y ); % transformed back 
    error_mae(i - 2002, 3) = 10 .^ median( abs( x - y ) );
    error_rmse(i - 2002, 3) = 10 .^ sqrt(mean( (x - y) .^ 2 ));
    c = corrcoef( x , y );
    error_r(i - 2002, 3) = c(2);
end

error_table_sp = table(years, ...
    error_mae(:,2), error_bias(:,2), error_r2(:,2), error_slope(:,2), error_r(:,2), ... % NEMO
    error_mae(:,1), error_bias(:,1), error_r2(:,1), error_slope(:,1), error_r(:,1), ... % ECCO
    error_mae(:,3), error_bias(:,3), error_r2(:,3), error_slope(:,3), error_r(:,3), ... % HYCOM
    'VariableNames', {'Year', ...
    'NEMO_MAE','NEMO_BIAS','NEMO_R2','NEMO_Slope','NEMO_r', ...
    'ECCO_MAE','ECCO_BIAS','ECCO_R2','ECCO_Slope','ECCO_r', ...
    'HYCOM_MAE','HYCOM_BIAS','HYCOM_R2','HYCOM_Slope','HYCOM_r'});

% error_table_sp{:, 2 : end} = round( error_table_sp{:,2:end}, 3 );
error_table_sp = sortrows( error_table_sp, 'Year', 'descend' );

%% Summer
for i = 2003 : 2022
    mld_argo = nan( 360, 180, 12, 'single' );
    for j = 6 : 8
            sub_folder1 = [PATH_ARGO, '\', '0.03sigmaTheta_ml_', num2str( i ), '.nc'];
            files = dir( [sub_folder1, '\', '*.nc'] );
            mld_n = ncread( [files(j).folder, '\', files(j).name], 'MLD' );
            mld_n(:, 1 : 90) = nan;
            mld_argo(:, :, j)  = mld_n;

            if j + 6 > 12
                k = j + 6 - 12;
            else
                k = j + 6;
            end
            mld_s = ncread( [files(k).folder, '\', files(k).name], 'MLD' );
            mld_s(:, 91 : 180) = nan;
            mld_argo(:, :, k)  = mld_s;
    end
    mld_argo = permute( mld_argo, [2, 1, 3] );

    d_ecco = load(  [PATH_ECCO, '\', 'mld_', num2str( i ), '.mat'] );
    d_nemo = load(  [PATH_NEMO, '\', 'mld_', num2str( i ), '.mat'] );
    d_hycom = load(  [PATH_HYCOM, '\', 'mld_', num2str( i ), '.mat'] );

    % mld_argo_array = [mld_argo_array; mld_argo(:)];
    % mld_nemo_array = [mld_nemo_array; d_nemo.mldG1deg(:)];
    % mld_ecco_array = [mld_ecco_array; d_ecco.mldG1deg(:)];
    % mld_hycom_array = [mld_hycom_array; d_hycom.mldG1deg(:)];

    mld_argo_array = [mld_argo(:)];
    mld_nemo_array = [d_nemo.mldG1deg(:)];
    mld_ecco_array = [d_ecco.mldG1deg(:)];
    mld_hycom_array = [d_hycom.mldG1deg(:)];

    % Argo vs ECCO
    idx  = ~isnan( mld_argo_array + mld_ecco_array );
    y = mld_argo_array(idx);
    x = mld_ecco_array(idx);

    x = log10( x ); 
    y = log10( y );
    N(i - 2002, 1) = length( y ); 
    mdl = fitlm( x, y );
    error_slope(i - 2002, 1) = mdl.Coefficients.Estimate(2);
    error_r2(i - 2002, 1) = mdl.Rsquared.Ordinary;
    error_bias(i - 2002, 1) = 10 .^ median( x - y ); % transformed back 
    error_mae(i - 2002, 1) = 10 .^ median( abs( x - y ) );
    error_rmse(i - 2002, 1) = 10 .^ sqrt(mean( (x - y) .^ 2 ));
    c = corrcoef( x , y );
    error_r(i - 2002, 1) = c(2);

    % Argo vs NEMO
    idx  = ~isnan( mld_argo_array + mld_nemo_array );
    y = mld_argo_array(idx);
    x = mld_nemo_array(idx);

    x = log10( x ); 
    y = log10( y );
    N(i - 2002, 2) = length( y ); 
    mdl = fitlm( x, y );
    error_slope(i - 2002, 2) = mdl.Coefficients.Estimate(2);
    error_r2(i - 2002, 2) = mdl.Rsquared.Ordinary;
    error_bias(i - 2002, 2) = 10 .^ median( x - y ); % transformed back 
    error_mae(i - 2002, 2) = 10 .^ median( abs( x - y ) );
    error_rmse(i - 2002, 2) = 10 .^ sqrt(mean( (x - y) .^ 2 ));
    c = corrcoef( x , y );
    error_r(i - 2002, 2) = c(2);

    % Argo vs HyCOM
    idx  = ~isnan( mld_argo_array + mld_hycom_array );
    y = mld_argo_array(idx);
    x = mld_hycom_array(idx);

    x = log10( x ); 
    y = log10( y );
    N(i - 2002, 3) = length( y ); 
    mdl = fitlm( x, y );
    error_slope(i - 2002, 3) = mdl.Coefficients.Estimate(2);
    error_r2(i - 2002, 3) = mdl.Rsquared.Ordinary;
    error_bias(i - 2002, 3) = 10 .^ median( x - y ); % transformed back 
    error_mae(i - 2002, 3) = 10 .^ median( abs( x - y ) );
    error_rmse(i - 2002, 3) = 10 .^ sqrt(mean( (x - y) .^ 2 ));
    c = corrcoef( x , y );
    error_r(i - 2002, 3) = c(2);
end

error_table_su = table(years, ...
    error_mae(:,2), error_bias(:,2), error_r2(:,2), error_slope(:,2), error_r(:,2), ... % NEMO
    error_mae(:,1), error_bias(:,1), error_r2(:,1), error_slope(:,1), error_r(:,1), ... % ECCO
    error_mae(:,3), error_bias(:,3), error_r2(:,3), error_slope(:,3), error_r(:,3), ... % HYCOM
    'VariableNames', {'Year', ...
    'NEMO_MAE','NEMO_BIAS','NEMO_R2','NEMO_Slope','NEMO_r', ...
    'ECCO_MAE','ECCO_BIAS','ECCO_R2','ECCO_Slope','ECCO_r', ...
    'HYCOM_MAE','HYCOM_BIAS','HYCOM_R2','HYCOM_Slope','HYCOM_r'});

% error_table_su{:, 2 : end} = round( error_table_su{:, 2 : end}, 3 );
error_table_su = sortrows( error_table_su, 'Year', 'descend' );

%% Autumn
for i = 2003 : 2022
    mld_argo = nan( 360, 180, 12, 'single' );
    for j = 9 : 11
            sub_folder1 = [PATH_ARGO, '\', '0.03sigmaTheta_ml_', num2str( i ), '.nc'];
            files = dir( [sub_folder1, '\', '*.nc'] );
            mld_n = ncread( [files(j).folder, '\', files(j).name], 'MLD' );
            mld_n(:, 1 : 90) = nan;
            mld_argo(:, :, j)  = mld_n;

            if j + 6 > 12
                k = j + 6 - 12;
            else
                k = j + 6;
            end
            mld_s = ncread( [files(k).folder, '\', files(k).name], 'MLD' );
            mld_s(:, 91 : 180) = nan;
            mld_argo(:, :, k)  = mld_s;
    end
    mld_argo = permute( mld_argo, [2, 1, 3] );

    d_ecco = load(  [PATH_ECCO, '\', 'mld_', num2str( i ), '.mat'] );
    d_nemo = load(  [PATH_NEMO, '\', 'mld_', num2str( i ), '.mat'] );
    d_hycom = load(  [PATH_HYCOM, '\', 'mld_', num2str( i ), '.mat'] );

    % mld_argo_array = [mld_argo_array; mld_argo(:)];
    % mld_nemo_array = [mld_nemo_array; d_nemo.mldG1deg(:)];
    % mld_ecco_array = [mld_ecco_array; d_ecco.mldG1deg(:)];
    % mld_hycom_array = [mld_hycom_array; d_hycom.mldG1deg(:)];

    mld_argo_array = [mld_argo(:)];
    mld_nemo_array = [d_nemo.mldG1deg(:)];
    mld_ecco_array = [d_ecco.mldG1deg(:)];
    mld_hycom_array = [d_hycom.mldG1deg(:)];

    % Argo vs ECCO
    idx  = ~isnan( mld_argo_array + mld_ecco_array );
    y = mld_argo_array(idx);
    x = mld_ecco_array(idx);

    x = log10( x ); 
    y = log10( y );
    N(i - 2002, 1) = length( y ); 
    mdl = fitlm( x, y );
    error_slope(i - 2002, 1) = mdl.Coefficients.Estimate(2);
    error_r2(i - 2002, 1) = mdl.Rsquared.Ordinary;
    error_bias(i - 2002, 1) = 10 .^ median( x - y ); % transformed back 
    error_mae(i - 2002, 1) = 10 .^ median( abs( x - y ) );
    error_rmse(i - 2002, 1) = 10 .^ sqrt(mean( (x - y) .^ 2 ));
    c = corrcoef( x , y );
    error_r(i - 2002, 1) = c(2);

    % Argo vs NEMO
    idx  = ~isnan( mld_argo_array + mld_nemo_array );
    y = mld_argo_array(idx);
    x = mld_nemo_array(idx);

    x = log10( x ); 
    y = log10( y );
    N(i - 2002, 2) = length( y ); 
    mdl = fitlm( x, y );
    error_slope(i - 2002, 2) = mdl.Coefficients.Estimate(2);
    error_r2(i - 2002, 2) = mdl.Rsquared.Ordinary;
    error_bias(i - 2002, 2) = 10 .^ median( x - y ); % transformed back 
    error_mae(i - 2002, 2) = 10 .^ median( abs( x - y ) );
    error_rmse(i - 2002, 2) = 10 .^ sqrt(mean( (x - y) .^ 2 ));
    c = corrcoef( x , y );
    error_r(i - 2002, 2) = c(2);

    % Argo vs HyCOM
    idx  = ~isnan( mld_argo_array + mld_hycom_array );
    y = mld_argo_array(idx);
    x = mld_hycom_array(idx);

    x = log10( x ); 
    y = log10( y );
    N(i - 2002, 3) = length( y ); 
    mdl = fitlm( x, y );
    error_slope(i - 2002, 3) = mdl.Coefficients.Estimate(2);
    error_r2(i - 2002, 3) = mdl.Rsquared.Ordinary;
    error_bias(i - 2002, 3) = 10 .^ median( x - y ); % transformed back 
    error_mae(i - 2002, 3) = 10 .^ median( abs( x - y ) );
    error_rmse(i - 2002, 3) = 10 .^ sqrt(mean( (x - y) .^ 2 ));
    c = corrcoef( x , y );
    error_r(i - 2002, 3) = c(2);
end

error_table_au = table(years, ...
    error_mae(:,2), error_bias(:,2), error_r2(:,2), error_slope(:,2), error_r(:,2), ... % NEMO
    error_mae(:,1), error_bias(:,1), error_r2(:,1), error_slope(:,1), error_r(:,1), ... % ECCO
    error_mae(:,3), error_bias(:,3), error_r2(:,3), error_slope(:,3), error_r(:,3), ... % HYCOM
    'VariableNames', {'Year', ...
    'NEMO_MAE','NEMO_BIAS','NEMO_R2','NEMO_Slope','NEMO_r', ...
    'ECCO_MAE','ECCO_BIAS','ECCO_R2','ECCO_Slope','ECCO_r', ...
    'HYCOM_MAE','HYCOM_BIAS','HYCOM_R2','HYCOM_Slope','HYCOM_r'});

% error_table_au{:, 2 : end} = round( error_table_au{:, 2 : end}, 3 );
error_table_au = sortrows( error_table_au,'Year','descend' );

%% Winter
for i = 2003 : 2022
    mld_argo = nan( 360, 180, 12, 'single' );
    for j = [1, 2, 12]
            sub_folder1 = [PATH_ARGO, '\', '0.03sigmaTheta_ml_', num2str( i ), '.nc'];
            files = dir( [sub_folder1, '\', '*.nc'] );
            mld_n = ncread( [files(j).folder, '\', files(j).name], 'MLD' );
            mld_n(:, 1 : 90) = nan;
            mld_argo(:, :, j)  = mld_n;

            if j + 6 > 12
                k = j + 6 - 12;
            else
                k = j + 6;
            end
            mld_s = ncread( [files(k).folder, '\', files(k).name], 'MLD' );
            mld_s(:, 91 : 180) = nan;
            mld_argo(:, :, k)  = mld_s;
    end
    mld_argo = permute( mld_argo, [2, 1, 3] );

    d_ecco = load(  [PATH_ECCO, '\', 'mld_', num2str( i ), '.mat'] );
    d_nemo = load(  [PATH_NEMO, '\', 'mld_', num2str( i ), '.mat'] );
    d_hycom = load(  [PATH_HYCOM, '\', 'mld_', num2str( i ), '.mat'] );

    % mld_argo_array = [mld_argo_array; mld_argo(:)];
    % mld_nemo_array = [mld_nemo_array; d_nemo.mldG1deg(:)];
    % mld_ecco_array = [mld_ecco_array; d_ecco.mldG1deg(:)];
    % mld_hycom_array = [mld_hycom_array; d_hycom.mldG1deg(:)];

    mld_argo_array = [mld_argo(:)];
    mld_nemo_array = [d_nemo.mldG1deg(:)];
    mld_ecco_array = [d_ecco.mldG1deg(:)];
    mld_hycom_array = [d_hycom.mldG1deg(:)];

    % Argo vs ECCO
    idx  = ~isnan( mld_argo_array + mld_ecco_array );
    y = mld_argo_array(idx);
    x = mld_ecco_array(idx);

    x = log10( x ); 
    y = log10( y );
    N(i - 2002, 1) = length( y ); 
    mdl = fitlm( x, y );
    error_slope(i - 2002, 1) = mdl.Coefficients.Estimate(2);
    error_r2(i - 2002, 1) = mdl.Rsquared.Ordinary;
    error_bias(i - 2002, 1) = 10 .^ median( x - y ); % transformed back 
    error_mae(i - 2002, 1) = 10 .^ median( abs( x - y ) );
    error_rmse(i - 2002, 1) = 10 .^ sqrt(mean( (x - y) .^ 2 ));
    c = corrcoef( x , y );
    error_r(i - 2002, 1) = c(2);

    % Argo vs NEMO
    idx  = ~isnan( mld_argo_array + mld_nemo_array );
    y = mld_argo_array(idx);
    x = mld_nemo_array(idx);

    x = log10( x ); 
    y = log10( y );
    N(i - 2002, 1) = length( y ); 
    mdl = fitlm( x, y );
    error_slope(i - 2002, 2) = mdl.Coefficients.Estimate(2);
    error_r2(i - 2002, 2) = mdl.Rsquared.Ordinary;
    error_bias(i - 2002, 2) = 10 .^ median( x - y ); % transformed back 
    error_mae(i - 2002, 2) = 10 .^ median( abs( x - y ) );
    error_rmse(i - 2002, 2) = 10 .^ sqrt(mean( (x - y) .^ 2 ));
    c = corrcoef( x , y );
    error_r(i - 2002, 2) = c(2);

    % Argo vs HyCOM
    idx  = ~isnan( mld_argo_array + mld_hycom_array );
    y = mld_argo_array(idx);
    x = mld_hycom_array(idx);

    x = log10( x ); 
    y = log10( y );
    N(i - 2002, 3) = length( y ); 
    mdl = fitlm( x, y );
    error_slope(i - 2002, 3) = mdl.Coefficients.Estimate(2);
    error_r2(i - 2002, 3) = mdl.Rsquared.Ordinary;
    error_bias(i - 2002, 3) = 10 .^ median( x - y ); % transformed back 
    error_mae(i - 2002, 3) = 10 .^ median( abs( x - y ) );
    error_rmse(i - 2002, 3) = 10 .^ sqrt(mean( (x - y) .^ 2 ));
    c = corrcoef( x , y );
    error_r(i - 2002, 3) = c(2);
end

error_table_wi = table(years, ...
    error_mae(:,2), error_bias(:,2), error_r2(:,2), error_slope(:,2), error_r(:,2), ... % NEMO
    error_mae(:,1), error_bias(:,1), error_r2(:,1), error_slope(:,1), error_r(:,1), ... % ECCO
    error_mae(:,3), error_bias(:,3), error_r2(:,3), error_slope(:,3), error_r(:,3), ... % HYCOM
    'VariableNames', {'Year', ...
    'NEMO_MAE','NEMO_BIAS','NEMO_R2','NEMO_Slope','NEMO_r', ...
    'ECCO_MAE','ECCO_BIAS','ECCO_R2','ECCO_Slope','ECCO_r', ...
    'HYCOM_MAE','HYCOM_BIAS','HYCOM_R2','HYCOM_Slope','HYCOM_r'});

% error_table_wi{:,2:end} = round( error_table_wi{:, 2 : end}, 2 );
error_table_wi = sortrows( error_table_wi, 'Year', 'descend' );

%%
filename = 'model_error_tables.xlsx';
writetable( error_table_sp, filename, 'Sheet','Spring' );
writetable( error_table_su, filename, 'Sheet','Summer' );
writetable( error_table_au, filename, 'Sheet','Autumn' );
writetable( error_table_wi, filename, 'Sheet','Winter' );

%%
% error_sp = varfun( @median, error_table_sp );
% error_sp.median_Year = [];
% 
% error_su = varfun( @median, error_table_su );
% error_su.median_Year = [];
% 
% error_au = varfun( @median, error_table_au );
% error_au.median_Year = [];
% 
% error_wi = varfun( @median, error_table_wi );
% error_wi.median_Year = [];

