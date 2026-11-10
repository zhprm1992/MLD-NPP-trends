% Used for: trends of the 0.125 MLD per latitude row and month, mld_trends_monthly_d125.mat (SI Fig. S2)
% Original: Manuscript\MLD\manuscript\nature\NCC\R1\code for revision\mld_d125\compute_mld_trends_monthly_d125.m (data paths updated to I:\Prj; computation unchanged)
% compute_mld_trends_monthly_d125.m
%
% Per-(latitude, calendar-month) Sen's-slope/Mann-Kendall MLD trend for
% the 0.125 kg/m^3 threshold, IDENTICAL method to the "mld-nemo" block of
% Code\npp_cbpm\fs6_heatmap_monthlyNPPtrends_ot.m (the script that
% produces mld_trends_pct_nemo / mld_trends_sig_nemo inside NPPtrends.mat,
% the actual input to the manuscript's real Figure 4 pipeline:
% NPPtrends.mat -> senMLD2_ot.m -> f4_plot_sensitivity_NPP_MLD.m). Same
% ktaub call (alpha=0.9999, wantplot=1), same dwtest autocorrelation
% check (computed but unused downstream, kept for fidelity), same %/yr
% normalization, same (12,nlat)->rot90->(nlat,12) shape convention -- only
% the input MLD source differs (mldLat_monthly_d125.mat instead of
% mldLat_monthly_nemo.mat).
%
% NPP is NOT recomputed here -- it doesn't depend on the MLD detection
% threshold, so npp_trends_pct_nemo / npp_trends_sig_nemo are reused
% unchanged from the existing Code\NPPtrends.mat by the downstream
% sensitivity script.
%
% Output: mld_trends_monthly_d125.mat (this folder) -- mld_trends_pct_nemo_d125,
% mld_trends_sig_nemo_d125 (both 120 lat x 12 month), lat.
%
% Haipeng Zhao
% 2026-07-26

clear; clc;
addpath('D:\Users\zhp19\OneDrive - Indiana State University\Research\MATLAB toolbox\ktaub');
DIR = fileparts(mfilename('fullpath'));

d = load(fullfile(DIR,'mldLat_monthly_d125.mat'));
idx = d.lat < 60 & d.lat > -60;
lat = d.lat(idx);
mldLat_nemo = d.mldLat(idx, :);

mld_trends_pct_nemo_d125 = nan( 12, length( lat ) );
mld_trends_nemo_d125     = nan( 12, length( lat ) );
mld_plotofslope_nemo_d125 = nan( 12, length( lat ) );
mld_trends_sig_nemo_d125  = nan( 12, length( lat ) );

p_ac_mld = nan(length(lat),1);
for i = 1 : length( lat )
    for j = 1 : 12
        mld_lat = mldLat_nemo(i, d.mm == j);
        if sum( isnan( mld_lat ) ) > 0
            continue
        end
        period_length = size( mld_lat, 2 );
        datain = [(1 : period_length)', mld_lat'];
        [~, ~, r] = regress( datain(:, 2), [(ones( length( datain(:, 1) ), 1 )), datain(:, 1)]  );
        [p_ac_mld(i), ~] = dwtest( r, datain(:, 1) );
        [~, ~, ~, sig, ~, ~, ~, sen, ~, plotofslope, ~] = ktaub(datain, 0.9999, 1);
        p_value = sig;
        if ~isnan(plotofslope)
            if plotofslope(1, 2) == 0
                mld_trends_pct_nemo_d125(j, i) = 0;
                mld_trends_nemo_d125(j, i) = 0;
                mld_plotofslope_nemo_d125(j, i) = 0;
                mld_trends_sig_nemo_d125(j, i) = p_value;
            else
                mld_trends_pct_nemo_d125(j, i) = single( (plotofslope(2, 2) - plotofslope(1, 2)) / plotofslope(1, 2) * 100 );
                mld_trends_nemo_d125(j, i) = sen;
                mld_plotofslope_nemo_d125(j, i) = plotofslope(1, 2);
                mld_trends_sig_nemo_d125(j, i) = p_value;
            end
        end
    end
    if mod(i,20)==0, fprintf('lat %d/%d done\n', i, length(lat)); end
end

mld_trends_pct_nemo_d125 = rot90( mld_trends_pct_nemo_d125 );
mld_trends_sig_nemo_d125 = rot90( mld_trends_sig_nemo_d125 );
mld_trends_nemo_d125     = rot90( mld_trends_nemo_d125 );
mld_plotofslope_nemo_d125 = rot90( mld_plotofslope_nemo_d125 );

save( fullfile(DIR,'mld_trends_monthly_d125.mat'), ...
    'mld_trends_pct_nemo_d125', 'mld_trends_sig_nemo_d125', ...
    'mld_trends_nemo_d125', 'mld_plotofslope_nemo_d125', 'lat' );
fprintf('Saved: %s\n', fullfile(DIR,'mld_trends_monthly_d125.mat'));
