% Used for: sensitivity with both thresholds, senMLD_v2_d125_vs_003.mat (Supplementary Fig. 2)
% Original: Manuscript\MLD\manuscript\nature\NCC\R1\code for revision\mld_d125\senMLD2_d125_vs_003_ot.m (data paths updated to I:\Prj; computation unchanged)
% senMLD2_d125_vs_003_ot.m
%
% NPP/MLD "sensitivity" (Sen's-slope %-trend ratio, significance-filtered
% at p<0.1 for BOTH trends, per 6 latitude bands x 4 seasons) computed
% IDENTICALLY to Code\senMLD2_ot.m, restricted to the NEMO product, run
% TWICE -- once with the manuscript's established 0.03 kg/m^3 MLD
% threshold (mld_trends_pct_nemo / mld_trends_sig_nemo, from the existing
% Code\NPPtrends.mat, unchanged), once with the Reviewer-1 0.125 kg/m^3
% sensitivity threshold (mld_trends_pct_nemo_d125 / mld_trends_sig_nemo_d125,
% from compute_mld_trends_monthly_d125.m). NPP trends are identical in
% both arms (npp_trends_pct_nemo / npp_trends_sig_nemo from NPPtrends.mat)
% since NPP does not depend on the MLD detection threshold -- only the
% MLD side of the ratio changes.
%
% Output: senMLD_v2_d125_vs_003.mat (this folder) -- senMLD_NEMO_003,
% senMLD_NEMO_d125 (each [6 lat-band x 4 season x 3 stat: median/std/n]),
% senPt_nemo_003, senPt_nemo_d125 (each {6x4} cell of raw [nppTrend,
% mldTrend] pairs, matching senMLD2_ot.m's senPt_nemo format exactly).
%
% Haipeng Zhao
% 2026-07-26

clear; clc;
DIR = fileparts(mfilename('fullpath'));

mmN = [3 4 5;
       6 7 8;
       9 10 11;
       12 1 2];  % [Spring, Summer, Autumn, Winter]

mmS = [9 10 11;
       12 1 2;
       3 4 5;
       6 7 8];

latBand = [60 40;
           40 20;
           20 0;
           0 -20;
           -20 -40;
           -40 -60];

senMLD_NEMO_003  = nan(6, 4, 3);
senMLD_NEMO_d125 = nan(6, 4, 3);

N = load(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), '4_trends_and_sensitivity', 'NPPtrends.mat'), ...
    'lat','npp_trends_pct_nemo','npp_trends_sig_nemo','mld_trends_pct_nemo','mld_trends_sig_nemo');
D125 = load(fullfile(DIR,'mld_trends_monthly_d125.mat'), 'lat','mld_trends_pct_nemo_d125','mld_trends_sig_nemo_d125');

assert(isequal(N.lat, D125.lat), '0.03 and 0.125 lat grids differ -- must match');
lat = flip( N.lat );

npp_trends_pct_nemo = N.npp_trends_pct_nemo;
npp_trends_sig_nemo = N.npp_trends_sig_nemo;

sources = { N.mld_trends_pct_nemo, N.mld_trends_sig_nemo, 'senMLD_NEMO_003', 'senPt_nemo_003'; ...
    D125.mld_trends_pct_nemo_d125, D125.mld_trends_sig_nemo_d125, 'senMLD_NEMO_d125', 'senPt_nemo_d125' };

for src = 1:2
    mld_trends_pct = sources{src,1};
    mld_trends_sig = sources{src,2};
    senMLD = nan(6,4,3);
    senPt = cell(6,4);

    for i = 1 : 6
        latIdx = lat < latBand(i, 1) & lat > latBand(i, 2);
        for j = 1 : 4
            c1_npp = nan( 120, 3 ); c1_npp_sig = nan( 120, 3 );
            c1_npp(lat > 0, :) = npp_trends_pct_nemo(lat > 0, mmN(j, :)); c1_npp_sig(lat > 0, :) = npp_trends_sig_nemo(lat > 0, mmN(j, :));
            c1_npp(lat < 0, :) = npp_trends_pct_nemo( lat < 0, mmS(j, :)); c1_npp_sig(lat < 0, :) = npp_trends_sig_nemo( lat < 0, mmS(j, :));
            c1_npp = c1_npp(latIdx, :); c1_npp_sig = c1_npp_sig(latIdx, :);
            c1_npp = c1_npp(:); c1_npp_sig = c1_npp_sig(:);
            idx_NPP_sig = ~isnan( c1_npp );
            c1_npp = c1_npp( idx_NPP_sig );
            c1_npp_sig = c1_npp_sig( idx_NPP_sig );

            c1_mld = nan( 120, 3 ); c1_mld_sig = nan( 120, 3 );
            c1_mld(lat > 0, :) = mld_trends_pct(lat > 0, mmN(j, :)); c1_mld_sig(lat > 0, :) = mld_trends_sig(lat > 0, mmN(j, :));
            c1_mld(lat < 0, :) = mld_trends_pct( lat < 0, mmS(j, :)); c1_mld_sig(lat < 0, :) = mld_trends_sig( lat < 0, mmS(j, :));
            c1_mld = c1_mld(latIdx, :); c1_mld_sig = c1_mld_sig(latIdx, :);
            c1_mld = c1_mld(:); c1_mld_sig = c1_mld_sig(:);
            c1_mld = c1_mld( idx_NPP_sig );
            c1_mld_sig = c1_mld_sig( idx_NPP_sig );

            sig_level = 0.1;
            idx_sig = (c1_npp_sig < sig_level) & (c1_mld_sig < sig_level);
            sen_NEMO = c1_npp(idx_sig) ./ c1_mld(idx_sig);
            if ~isempty( sen_NEMO )
                senMLD(i, j, 1) = nanmedian( sen_NEMO(:, 1) );
                senMLD(i, j, 2) = nanstd( sen_NEMO(:, 1) );
                senMLD(i, j, 3) = sum( ~isnan( sen_NEMO(:, 1) ) );
            end
            c1_mld(~idx_sig) = 0;
            senPt{i, j} = [c1_npp, c1_mld];
        end
    end

    assignin('base', sources{src,3}, senMLD);
    assignin('base', sources{src,4}, senPt);
    fprintf('%s: %d/%d lat-band x season cells have n>=1\n', sources{src,3}, sum(~isnan(senMLD(:,:,1)),'all'), 24);
end

save( fullfile(DIR,'senMLD_v2_d125_vs_003.mat'), ...
    'senMLD_NEMO_003', 'senMLD_NEMO_d125', 'senPt_nemo_003', 'senPt_nemo_d125' );
fprintf('Saved: %s\n', fullfile(DIR,'senMLD_v2_d125_vs_003.mat'));
