% Used for: per-grid-cell monthly trends of NPP-CbPM and MLD-NEMO, grid_monthly_trends_mld_npp_nemo.mat (Extended Data Fig. 9)
% Original: Manuscript\MLD\manuscript\nature\NCC\R1\code for revision\mld_npp_relationship_biome\compute_grid_monthly_trends.m (data paths updated to I:\Prj; computation unchanged)
% compute_grid_monthly_trends.m
%
% Stage 1 of the biome-scale sensitivity analysis (biome-scale analog of
% Fig. 4 / senMLD2_ot.m, at grid-cell resolution instead of lat-row
% resolution).
%
% For every 1x1 deg ocean grid cell, and separately for each of the 12
% calendar months, compute a 20-year (2003-2022) Sen's slope trend (%,
% ktaub) for NPP and for MLD, with significance -- mirroring
% Code\npp_cafe\heatmap_monthlyNPPtrends_ot.m's per-lat-row-per-month
% method exactly (same ktaub call, same trend_pct formula), just applied
% to every (lat,lon) grid cell instead of a lat-row zonal mean.
%
% Data sources (matches the CbPM/NEMO combination used throughout this R1
% revision, NOT the original Fig 4 script's CAFE + stale-MLD-path
% combination -- confirmed with user 2026-07-30):
%   NPP (CbPM, zno3mMU): I:\Prj\npp_cbpm_1M_1dg_zno3mMU\grid_npp_nemo\npp_YYYY.mat  (nppG1dg)
%   MLD (NEMO, corrected one_d03 source): I:\Prj\MLD_1M_NEMO\one_d03\mld\mld_YYYY.mat  (mldG1deg)
%
% A cell/month is only scored if all 20 years are finite for that
% variable (same completeness requirement as the lat-row version).
%
% Output: grid_monthly_trends_mld_npp_nemo.mat (this folder) --
%   npp_trend_pct, npp_trend_p, mld_trend_pct, mld_trend_p : (180,360,12)
%   lat (180,1), lon (360,1)
%
% Haipeng Zhao
% 2026-07-30

clear; clc;
addpath('D:\Users\zhp19\OneDrive - Indiana State University\Research\MATLAB toolbox\ktaub');

NPP_DIR = 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\grid_npp_nemo';
MLD_DIR = 'I:\Prj\MLD_1M_NEMO\one_d03\mld';
YEARS = 2003:2022;
NYR = numel(YEARS);
NLAT = 180; NLON = 360; NMO = 12;

ref = load(fullfile(NPP_DIR, sprintf('npp_%d.mat', YEARS(1))), 'lat','lon');
lat = double(ref.lat(:)); lon = double(ref.lon(:));

fprintf('Loading %d years of NPP and MLD...\n', NYR); t0 = tic;
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
save(fullfile(OUT_DIR, 'grid_monthly_trends_mld_npp_nemo.mat'), ...
    'npp_trend_pct','npp_trend_p','mld_trend_pct','mld_trend_p','lat','lon','-v7.3');
fprintf('Saved: %s\n', fullfile(OUT_DIR, 'grid_monthly_trends_mld_npp_nemo.mat'));
