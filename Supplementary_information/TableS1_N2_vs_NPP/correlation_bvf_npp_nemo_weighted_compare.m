% Used for: SI Table S1 values (the "new" depth-weighted columns)
% Original: Manuscript\MLD\manuscript\nature\NCC\R1\code for revision\correlation_bvf_npp\correlation_bvf_npp_nemo_weighted_compare.m (data paths updated to I:\Prj; computation unchanged)
% correlation_bvf_npp_nemo_weighted_compare.m
%
% Reruns correlation_bvf_npp_nemo_run.m's exact analysis (the script
% behind Table S5/S6 in the R1 response letter) with the new
% depth-weighted 0-200m BVF (bvfLat_season_weighted.mat, built via
% prj_grid_bvf200_weighted_ot.m -> biomes_1M_bvf200_weighted_ot.m ->
% latSeason_bvf200_weighted_ot.m, fixing the np.trapezoid/numpy-1.24 bug
% and using a depth-weighted trapezoidal average instead of the old
% unweighted level-mean), side by side with the original
% bvfLat_season.mat, to see whether the depth-weighting revision changes
% the reported NPP-BVF correlation/trend results.
%
% Same 6 lat bands x 4 seasons, same Sen's-slope trend method, same
% Pearson correlation -- only the BVF source data differs between the two
% columns.
%
% Haipeng Zhao
% 2026-07-29

clear
clc
DIR = fileparts(mfilename('fullpath'));
addpath('D:\Users\zhp19\OneDrive - Indiana State University\Research\MATLAB toolbox\ktaub');

latBand = [60 40;
           40 20;
           20 0;
           0 -20;
           -20 -40;
           -40 -60];

seasons = {'sp', 'su', 'au', 'wi'};
seasonNames = {'Spring', 'Summer', 'Autumn', 'Winter'};

d_bvf_old = load( 'I:\Prj\MLD_1M_NEMO\bvfLat_season.mat' );
d_bvf_new = load( fullfile(DIR, 'bvfLat_season_weighted.mat') );
d_npp = load( 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\nppLat_season_nemo.mat' );

lat = d_bvf_old.lat(:);
L = 20;

nBand = size( latBand, 1 );
nSeason = length( seasons );

sources = {d_bvf_old, 'old'; d_bvf_new, 'new'};
results = struct();

for src = 1:2
    d_bvf = sources{src,1}; tag = sources{src,2};
    npp_trend = nan(nBand,nSeason); bvf_trend = nan(nBand,nSeason);
    npp_trend_sig = nan(nBand,nSeason); bvf_trend_sig = nan(nBand,nSeason);
    r_corr = nan(nBand,nSeason); p_corr = nan(nBand,nSeason);

    for j = 1 : nSeason
        bvfLat = d_bvf.( ['bvf', seasons{j}, 'Lat'] );
        nppLat = d_npp.( ['npp', seasons{j}, 'Lat'] );

        for i = 1 : nBand
            latIdx = lat < latBand(i, 1) & lat > latBand(i, 2);

            bvf_ts = nanmean( bvfLat(latIdx, :), 1 )';
            npp_ts = nanmean( nppLat(latIdx, :), 1 )';

            datain = [(1:L)', npp_ts];
            [~,~,~,sig,~,~,~,~,~,plotofslope,~] = ktaub(datain, 0.9999, 1);
            if ~isnan(plotofslope)
                npp_trend(i,j) = single((plotofslope(2,2)-plotofslope(1,2))/plotofslope(1,2)*100);
            end
            npp_trend_sig(i,j) = sig;

            datain = [(1:L)', bvf_ts];
            [~,~,~,sig,~,~,~,~,~,plotofslope,~] = ktaub(datain, 0.9999, 1);
            if ~isnan(plotofslope)
                bvf_trend(i,j) = single((plotofslope(2,2)-plotofslope(1,2))/abs(plotofslope(1,2))*100);
            end
            bvf_trend_sig(i,j) = sig;

            [R,P] = corrcoef(npp_ts, bvf_ts, 'rows','complete');
            r_corr(i,j) = R(1,2); p_corr(i,j) = P(1,2);
        end
    end
    results.(tag).npp_trend = npp_trend;
    results.(tag).bvf_trend = bvf_trend;
    results.(tag).npp_trend_sig = npp_trend_sig;
    results.(tag).bvf_trend_sig = bvf_trend_sig;
    results.(tag).r_corr = r_corr;
    results.(tag).p_corr = p_corr;
end

fprintf('%-12s %-8s | %10s %10s | %10s %10s | %10s %10s\n', 'Band','Season', ...
    'BVFtr_old','BVFtr_new','r_old','r_new','p_old','p_new');
for j = 1:nSeason
    for i = 1:nBand
        sign_flip = sign(results.old.r_corr(i,j)) ~= sign(results.new.r_corr(i,j));
        flag = ''; if sign_flip, flag = '  <-- r sign flip'; end
        fprintf('%3d to %3d  %-8s | %10.3f %10.3f | %+10.3f %+10.3f | %10.4f %10.4f%s\n', ...
            latBand(i,1), latBand(i,2), seasonNames{j}, ...
            results.old.bvf_trend(i,j), results.new.bvf_trend(i,j), ...
            results.old.r_corr(i,j), results.new.r_corr(i,j), ...
            results.old.p_corr(i,j), results.new.p_corr(i,j), flag);
    end
end

out_mat = fullfile(DIR, 'corr_bvf_npp_nemo_weighted_compare.mat');
save(out_mat, 'results', 'latBand', 'seasons', 'seasonNames');
fprintf('\nSaved: %s\n', out_mat);
