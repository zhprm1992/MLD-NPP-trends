% Used for: per-grid-cell monthly trends of ensemble NPP and MLD-NEMO, grid_monthly_trends_mld_npp_ensemble_nemo.mat (Supplementary Fig. 5)
% Original: Manuscript\MLD\manuscript\nature\NCC\R1\code for revision\mld_npp_relationship_biome\compute_grid_monthly_trends_ensemble.m (data paths updated to I:\Prj; computation unchanged)
% compute_grid_monthly_trends_ensemble.m
%
% Ensemble-NPP counterpart to compute_grid_monthly_trends.m. Identical
% method (per-cell, per-calendar-month, 20-yr Theil-Sen/ktaub trend with
% significance), only the NPP source is swapped from CbPM to the
% VGPM+CAFE+CbPM ensemble mean. MLD source is unchanged (same MLD-NEMO
% series underlies both the CbPM and ensemble comparisons).
%
% Data sources:
%   NPP (VGPM+CAFE+CbPM ensemble): I:\Prj\npp_ensemble_1M_1dg\grid_npp_nemo_1M_1dg_mat\npp_YYYY.mat  (nppG1dg)
%   MLD (NEMO, corrected one_d03 source): I:\Prj\MLD_1M_NEMO\one_d03\mld\mld_YYYY.mat  (mldG1deg)
%
% A cell/month is only scored if all 20 years are finite for that
% variable (same completeness requirement as the CbPM version).
%
% Output: grid_monthly_trends_mld_npp_ensemble_nemo.mat (this folder) --
%   npp_trend_pct, npp_trend_p, mld_trend_pct, mld_trend_p : (180,360,12)
%   lat (180,1), lon (360,1)
%
% Haipeng Zhao
% 2026-09-12

clear; clc;
addpath('D:\Users\zhp19\OneDrive - Indiana State University\Research\MATLAB toolbox\ktaub');

NPP_DIR = 'I:\Prj\npp_ensemble_1M_1dg\grid_npp_nemo_1M_1dg_mat';
MLD_DIR = 'I:\Prj\MLD_1M_NEMO\one_d03\mld';
YEARS = 2003:2022;
NYR = numel(YEARS);
NLAT = 180; NLON = 360; NMO = 12;

ref = load(fullfile(NPP_DIR, sprintf('npp_%d.mat', YEARS(1))), 'lat','lon');
lat = double(ref.lat(:)); lon = double(ref.lon(:));

fprintf('Loading %d years of ensemble NPP and MLD...\n', NYR); t0 = tic;
npp_mo = nan(NLAT,NLON,NMO,NYR,'single');
mld_mo = nan(NLAT,NLON,NMO,NYR,'single');
for y = 1:NYR
    d = load(fullfile(NPP_DIR, sprintf('npp_%d.mat', YEARS(y))), 'nppG1dg');
    npp_mo(:,:,:,y) = d.nppG1dg;
    d = load(fullfile(MLD_DIR, sprintf('mld_%d.mat', YEARS(y))), 'mldG1deg');
    mld_mo(:,:,:,y) = single(d.mldG1deg);
end
fprintf('Loaded in %.0fs\n', toc(t0));

npp_trend_pct = nan(NLAT,NLON,NMO,'single');
npp_trend_p   = nan(NLAT,NLON,NMO,'single');
mld_trend_pct = nan(NLAT,NLON,NMO,'single');
mld_trend_p   = nan(NLAT,NLON,NMO,'single');

x = (1:NYR)';
t0 = tic;
for m = 1:NMO
    fprintf('Month %d/12 ...\n', m); tm = tic;
    for i = 1:NLAT
        for j = 1:NLON
            y_npp = squeeze(npp_mo(i,j,m,:));
            if sum(isnan(y_npp)) == 0
                datain = [x, y_npp];
                [~,~,~,sig,~,~,~,~,~,plotofslope,~] = ktaub(datain, 0.9999, 1);
                if ~isnan(plotofslope)
                    if plotofslope(1,2) == 0
                        npp_trend_pct(i,j,m) = 0;
                    else
                        npp_trend_pct(i,j,m) = single((plotofslope(2,2)-plotofslope(1,2))/plotofslope(1,2)*100);
                    end
                    npp_trend_p(i,j,m) = sig;
                end
            end

            y_mld = squeeze(mld_mo(i,j,m,:));
            if sum(isnan(y_mld)) == 0
                datain = [x, y_mld];
                [~,~,~,sig,~,~,~,~,~,plotofslope,~] = ktaub(datain, 0.9999, 1);
                if ~isnan(plotofslope)
                    if plotofslope(1,2) == 0
                        mld_trend_pct(i,j,m) = 0;
                    else
                        mld_trend_pct(i,j,m) = single((plotofslope(2,2)-plotofslope(1,2))/plotofslope(1,2)*100);
                    end
                    mld_trend_p(i,j,m) = sig;
                end
            end
        end
    end
    fprintf('  month %d done in %.0fs\n', m, toc(tm));
end
fprintf('All months done in %.0fs\n', toc(t0));

OUT_DIR = fileparts(mfilename('fullpath'));
save(fullfile(OUT_DIR, 'grid_monthly_trends_mld_npp_ensemble_nemo.mat'), ...
    'npp_trend_pct','npp_trend_p','mld_trend_pct','mld_trend_p','lat','lon','-v7.3');
fprintf('Saved: %s\n', fullfile(OUT_DIR, 'grid_monthly_trends_mld_npp_ensemble_nemo.mat'));
