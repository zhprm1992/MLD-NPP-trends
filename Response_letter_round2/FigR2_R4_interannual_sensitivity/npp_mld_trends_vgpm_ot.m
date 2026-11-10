% Used for: VGPM trends per latitude row and month, VGPM_NPPtrends.mat (Figs. R2, R4)
% Original: Manuscript\MLD\manuscript\nature\NCC\R1\code for revision\npp_ensemble\npp_mld_trends_vgpm_ot.m (data paths updated to I:\Prj; computation unchanged)
% Compute Mann-Kendall trends (% yr-1) of NPP (VGPM) and MLD (NEMO) by
% latitude and month, for the NPP-vs-MLD sensitivity figure.
% Adapted from code for revision/npp_cafe/heatmap_monthlyNPPtrends_ot.m
% (computation only, NEMO-only, no heatmap plotting).
% Haipeng Zhao
% 7/12/2026

clear
clc

d = load( 'I:\Prj\npp_vgpm_1M_1dg\nppLat_monthly_vgpm.mat' );
idx = d.lat < 60 & d.lat > -60;
lat = d.lat(idx);
nppLat = d.nppLat(idx, :);
mm_npp = d.mm;

% Align to the MLD-NEMO period (2003-2022); VGPM's lat-profile has an
% extra year (2003-2023) with no MLD-NEMO counterpart yet.
nYearsMLD = 20;
nppLat = nppLat(:, 1 : nYearsMLD * 12);
mm_npp = mm_npp(1 : nYearsMLD * 12);

d2 = load( 'I:\Prj\MLD_1M_NEMO\mldLat_monthly_nemo.mat' );
idx2 = d2.lat < 60 & d2.lat > -60;
mldLat = d2.mldLat(idx2, :);
mm_mld = d2.mm;

npp_trends_pct_nemo = nan( 12, length( lat ) );
npp_trends_nemo = nan( 12, length( lat ) );
npp_plotofslope_nemo = nan( 12, length( lat ) );
npp_trends_sig_nemo = nan( 12, length( lat ) );

mld_trends_pct_nemo = nan( 12, length( lat ) );
mld_trends_nemo = nan( 12, length( lat ) );
mld_plotofslope_nemo = nan( 12, length( lat ) );
mld_trends_sig_nemo = nan( 12, length( lat ) );

%% NPP trends
for i = 1 : length( lat )
    for j = 1 : 12
        npp_lat = nppLat(i, mm_npp == j);
        if sum( isnan( npp_lat ) ) > 0
            continue
        end
        period_length = size( npp_lat, 2 );
        datain = [(1 : period_length)', npp_lat'];
        [~, ~, r] = regress( datain(:, 2), [(ones( length( datain(:, 1) ), 1 )), datain(:, 1)]  );
        dwtest( r, datain(:, 1) );
        [~, ~, ~, sig, ~, ~, ~, sen, ~, plotofslope, ~] = ktaub(datain, 0.9999, 1);
        p_value = sig;
        if ~isnan(plotofslope)
            if plotofslope(1, 2) == 0
                npp_trends_pct_nemo(j, i) = 0;
                npp_trends_nemo(j, i) = 0;
                npp_plotofslope_nemo(j, i) = 0;
                npp_trends_sig_nemo(j, i) = p_value;
            else
                npp_trends_pct_nemo(j, i) = single( (plotofslope(2, 2) - plotofslope(1, 2)) / plotofslope(1, 2) * 100 );
                npp_trends_nemo(j, i) = sen;
                npp_plotofslope_nemo(j, i) = plotofslope(1, 2);
                npp_trends_sig_nemo(j, i) = p_value;
            end
        end
    end
end

%% MLD trends
for i = 1 : length( lat )
    for j = 1 : 12
        mld_lat = mldLat(i, mm_mld == j);
        if sum( isnan( mld_lat ) ) > 0
            continue
        end
        period_length = size( mld_lat, 2 );
        datain = [(1 : period_length)', mld_lat'];
        [~, ~, r] = regress( datain(:, 2), [(ones( length( datain(:, 1) ), 1 )), datain(:, 1)]  );
        dwtest( r, datain(:, 1) );
        [~, ~, ~, sig, ~, ~, ~, sen, ~, plotofslope, ~] = ktaub(datain, 0.9999, 1);
        p_value = sig;
        if ~isnan(plotofslope)
            if plotofslope(1, 2) == 0
                mld_trends_pct_nemo(j, i) = 0;
                mld_trends_nemo(j, i) = 0;
                mld_plotofslope_nemo(j, i) = 0;
                mld_trends_sig_nemo(j, i) = p_value;
            else
                mld_trends_pct_nemo(j, i) = single( (plotofslope(2, 2) - plotofslope(1, 2)) / plotofslope(1, 2) * 100 );
                mld_trends_nemo(j, i) = sen;
                mld_plotofslope_nemo(j, i) = plotofslope(1, 2);
                mld_trends_sig_nemo(j, i) = p_value;
            end
        end
    end
end

% Rotate to [lat x month], matching the orientation expected by
% senMLD_vgpm_ot.m (mirrors the rot90() step applied before the heatmap
% plot / save in the CAFE original)
npp_trends_pct_nemo = rot90( npp_trends_pct_nemo );
npp_trends_sig_nemo = rot90( npp_trends_sig_nemo );
mld_trends_pct_nemo = rot90( mld_trends_pct_nemo );
mld_trends_sig_nemo = rot90( mld_trends_sig_nemo );

save( 'I:\Prj\npp_vgpm_1M_1dg\VGPM_NPPtrends.mat', 'npp_trends_pct_nemo', 'mld_trends_pct_nemo',...
      'npp_trends_sig_nemo', 'mld_trends_sig_nemo', 'lat' )
disp( 'Save successfully: VGPM_NPPtrends' )
