% Used for: Extended Data Fig. 6 radar chart (reads error_metrics_season.mat, the four seasonal tables of scatterMLD_ARGOvsModels_point_season.m; needs radar-chart toolbox, see README)
% Original: Code\spiderplot_v3.m (data paths updated to I:\Prj; computation unchanged)
% Create star  plot of error metrics (mae, bias, rmse, r-square, r), 
% Haipeng Zhao
% 1/15/2025
% 2/17/2025

clear
clc

fig = figure;
fig.Position = [100 100 900 700];  % [left bottom width height] in pixels
t = tiledlayout( 2, 2 );
t.Padding  = "compact" ;
t.TileSpacing = "loose";

load error_metrics_season.mat

%% Spring
nexttile
% MAE, Bias, r, r-square, slope
mae_score = 1 - abs( error_sp.median_MAE_NEMO - 1 );
bias_score = 1 - abs ( error_sp.median_Bias_NEMO  - 1 );
slope_score = 1 - abs ( error_sp.median_Slope_NEMO  - 1 );
nemo_metrics = [mae_score, bias_score, error_sp.median_r_NEMO, error_sp.median_R2_NEMO, slope_score];

mae_score = 1 - abs( error_sp.median_MAE_ECCO - 1 );
bias_score = 1 - abs ( error_sp.median_Bias_ECCO  - 1 );
slope_score = 1 - abs ( error_sp.median_Slope_ECCO  - 1 );
ecco_metrics = [mae_score, bias_score, error_sp.median_r_ECCO, error_sp.median_R2_ECCO, slope_score];

mae_score = 1 - abs( error_sp.median_MAE_HYCOM - 1 );
bias_score = 1 - abs ( error_sp.median_Bias_HYCOM  - 1 );
slope_score = 1 - abs ( error_sp.median_Slope_HYCOM  - 1 );
hycom_metrics = [mae_score, bias_score, error_sp.median_r_HYCOM, error_sp.median_R2_HYCOM, slope_score];

P = [nemo_metrics; ecco_metrics; hycom_metrics];
radar_chart( P )

%% Summer
nexttile
mae_score = 1 - abs( error_su.median_MAE_NEMO - 1 );
bias_score = 1 - abs ( error_su.median_Bias_NEMO  - 1 );
slope_score = 1 - abs ( error_su.median_Slope_NEMO  - 1 );
nemo_metrics = [mae_score, bias_score, error_su.median_r_NEMO, error_su.median_R2_NEMO, slope_score];

mae_score = 1 - abs( error_su.median_MAE_ECCO - 1 );
bias_score = 1 - abs ( error_su.median_Bias_ECCO  - 1 );
slope_score = 1 - abs ( error_su.median_Slope_ECCO  - 1 );
ecco_metrics = [mae_score, bias_score, error_su.median_r_ECCO, error_su.median_R2_ECCO, slope_score];

mae_score = 1 - abs( error_su.median_MAE_HYCOM - 1 );
bias_score = 1 - abs ( error_su.median_Bias_HYCOM  - 1 );
slope_score = 1 - abs ( error_su.median_Slope_HYCOM  - 1 );
hycom_metrics = [mae_score, bias_score, error_su.median_r_HYCOM, error_su.median_R2_HYCOM, slope_score];

P = [nemo_metrics; ecco_metrics; hycom_metrics];
radar_chart( P )


%% Autumn
nexttile
mae_score = 1 - abs( error_au.median_MAE_NEMO - 1 );
bias_score = 1 - abs ( error_au.median_Bias_NEMO  - 1 );
slope_score = 1 - abs ( error_au.median_Slope_NEMO  - 1 );
nemo_metrics = [mae_score, bias_score, error_au.median_r_NEMO, error_au.median_R2_NEMO, slope_score];

mae_score = 1 - abs( error_au.median_MAE_ECCO - 1 );
bias_score = 1 - abs ( error_au.median_Bias_ECCO  - 1 );
slope_score = 1 - abs ( error_au.median_Slope_ECCO  - 1 );
ecco_metrics = [mae_score, bias_score, error_au.median_r_ECCO, error_au.median_R2_ECCO, slope_score];

mae_score = 1 - abs( error_au.median_MAE_HYCOM - 1 );
bias_score = 1 - abs ( error_au.median_Bias_HYCOM  - 1 );
slope_score = 1 - abs ( error_au.median_Slope_HYCOM  - 1 );
hycom_metrics = [mae_score, bias_score, error_au.median_r_HYCOM, error_au.median_R2_HYCOM, slope_score];

P = [nemo_metrics; ecco_metrics; hycom_metrics];
radar_chart( P )

%% Winter
nexttile
mae_score = 1 - abs( error_wi.median_MAE_NEMO - 1 );
bias_score = 1 - abs ( error_wi.median_Bias_NEMO  - 1 );
slope_score = 1 - abs ( error_wi.median_Slope_NEMO  - 1 );
nemo_metrics = [mae_score, bias_score, error_wi.median_r_NEMO, error_wi.median_R2_NEMO, slope_score];

mae_score = 1 - abs( error_wi.median_MAE_ECCO - 1 );
bias_score = 1 - abs ( error_wi.median_Bias_ECCO  - 1 );
slope_score = 1 - abs ( error_wi.median_Slope_ECCO  - 1 );
ecco_metrics = [mae_score, bias_score, error_wi.median_r_ECCO, error_wi.median_R2_ECCO, slope_score];

mae_score = 1 - abs( error_wi.median_MAE_HYCOM - 1 );
bias_score = 1 - abs ( error_wi.median_Bias_HYCOM  - 1 );
slope_score = 1 - abs ( error_wi.median_Slope_HYCOM  - 1 );
hycom_metrics = [mae_score, bias_score, error_wi.median_r_HYCOM, error_wi.median_R2_HYCOM, slope_score];

P = [nemo_metrics; ecco_metrics; hycom_metrics];
radar_chart( P )

