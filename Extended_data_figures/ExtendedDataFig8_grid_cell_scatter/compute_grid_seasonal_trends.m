% Used for: per-grid-cell seasonal trends of NPP-CbPM and MLD-NEMO, grid_trends_mld_npp.mat (Extended Data Fig. 8)
% Original: Manuscript\MLD\manuscript\nature\NCC\R1\code for revision\mld_npp_relationship_grid\compute_grid_seasonal_trends.m (data paths updated to I:\Prj; computation unchanged)
% compute_grid_seasonal_trends.m
%
% Grid-cell (1x1 deg) seasonal trend analysis of MLD and CbPM NPP, replacing
% the latitudinal-bin approach in ..\correlation_bvf_npp and ..\npp_cafe with
% a true 2-D grid-cell analysis, analogous to the Cael et al. (2026) Fig. 2
% grid-cell design in I:\Prj3\Cael_2026\cael_replication_and_correction.m.
%
% For each 1x1 deg grid cell, each of the 3 MLD/NPP products (NEMO, ECCO,
% HYCOM-gat), and each of the 4 hemisphere-aware seasons:
%   1. Build the 20 seasonal-mean values (one per year, 2003-2022) from the
%      monthly grids.
%   2. Fit an OLS linear trend vs. centered year index (t = -9.5:9.5), giving
%      slope b (+-SE) and intercept a (+-SE, = mean state since t is centered).
%   3. Convert to a relative trend: trend_pct = b/a*100 (%/yr), with the
%      CORRECTED uncertainty propagation (matching the fix identified in
%      Cael's own analysis, not the original paper's buggy formula):
%        trend_pct_se = |trend_pct| * sqrt((SEb/b)^2 + (SEa/a)^2)
%   A cell is only scored if all 20 seasonal-mean years are finite (same
%   completeness requirement Cael used for the open-ocean mask).
%
% Season definition (hemisphere-shifted, calendar-year pooled -- matches the
% convention already used in ..\npp_cafe\senMLD2_ot.m and Prj4\Code\Chl_LOV\
% chl3d_1S_ot.m: no Dec->next-year wraparound, each row uses its own
% hemisphere's month grouping from the SAME calendar-year file):
%   mmN = [Spring MAM; Summer JJA; Autumn SON; Winter DJF]  (lat >= 0)
%   mmS = [Spring SON; Summer DJF; Autumn MAM; Winter JJA]  (lat <  0)
%
% Inputs (read-only, on I: drive -- not modified):
%   NPP  (CbPM, zno3mMU variant): I:\Prj\npp_cbpm_1M_1dg_zno3mMU\grid_npp_{nemo,ecco,hycom}\npp_YYYY.mat  (nppG1dg)
%   MLD                         : I:\Prj\MLD_1M_NEMO\one_d03\mld\mld_YYYY.mat                             (mldG1deg)
%                                  I:\Prj\MLD_ECCO\one\mld_YYYY.mat                                       (mldG1deg)
%                                  I:\Prj\MLD_HYCOM\one_gat\mld\mld_YYYY.mat                              (mldG1deg)
%
% CORRECTED 2026-07-25: NEMO MLD source was I:\Prj\MLD_1M_NEMO\one\, a
% stale (March 2025) product that does NOT match what CbPM's own NPP
% generator (Code\npp_cbpm\CbPM_mld.m) or the manuscript's own validation
% script (scatterMLD_ARGOvsMODEL_trends_sp.m) actually use -- both source
% MLD from I:\Prj\MLD_1M_NEMO\one_d03\mld\ (May 2025). The two products
% differ substantially (median 28.85 vs 38.58 m, ~1/777600 cells match
% exactly). Since CbPM's NPP is itself computed using MLD (photoacclimation/
% depth-integration), pairing it with a different MLD field than the one it
% was actually computed against was a real provenance bug. User confirmed
% one_d03 is correct; one\ must not be used anywhere in this analysis.
% Impact of the fix on the final biome-level results: 54/55 (98%) sign
% agreement with the pre-fix values, r=0.96, median change ~0 -- the
% correction was necessary but did not change the qualitative conclusion.
%
% Output: grid_trends_mld_npp.mat (this folder) -- lat, lon, season_names,
% product_names, and NPP/MLD trend_pct / trend_se / mean_state / n_valid
% grids, size (180,360,4,3).
%
% Haipeng Zhao

clear; close all; clc;
OUT_DIR = fullfile(fileparts(mfilename('fullpath')));

YEARS   = 2003:2022;
NYR     = numel(YEARS);
NLAT    = 180; NLON = 360; NMO = 12;

product_names = {'nemo','ecco','hycom'};
npp_dir = struct( ...
    'nemo',  'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\grid_npp_nemo', ...
    'ecco',  'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\grid_npp_ecco', ...
    'hycom', 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\grid_npp_hycom');
mld_dir = struct( ...
    'nemo',  'I:\Prj\MLD_1M_NEMO\one_d03\mld', ...   % CORRECTED 2026-07-25, was 'I:\Prj\MLD_1M_NEMO\one'
    'ecco',  'I:\Prj\MLD_ECCO\one', ...
    'hycom', 'I:\Prj\MLD_HYCOM\one_gat\mld');

season_names = {'Spring','Summer','Autumn','Winter'};
mmN = [3 4 5; 6 7 8; 9 10 11; 12 1 2];   % lat >= 0
mmS = [9 10 11; 12 1 2; 3 4 5; 6 7 8];   % lat <  0

% Reference grid
ref = load(fullfile(npp_dir.nemo, sprintf('npp_%d.mat', YEARS(1))), 'lat','lon');
lat = double(ref.lat(:));   % (180,1), -89.5 .. 89.5
lon = double(ref.lon(:));   % (360,1),   0.5 .. 359.5
assert(numel(lat)==NLAT && numel(lon)==NLON);
idxN = lat >= 0;
idxS = lat <  0;

t_vec = (-9.5:9.5)';                 % centered year index, NYR=20
Sxx   = sum(t_vec.^2);
dof   = NYR - 2;

NPP.trend_pct  = nan(NLAT,NLON,4,3,'single');
NPP.trend_se   = nan(NLAT,NLON,4,3,'single');
NPP.mean_state = nan(NLAT,NLON,4,3,'single');
NPP.n_valid    = zeros(NLAT,NLON,4,3,'uint8');
MLD = NPP;

for p = 1:3
    pname = product_names{p};
    fprintf('\n=== Product %d/3: %s ===\n', p, upper(pname));

    fprintf('  Loading %d years of NPP ...\n', NYR); t0 = tic;
    npp_mo = nan(NLAT,NLON,NMO,NYR,'single');
    for y = 1:NYR
        f = fullfile(npp_dir.(pname), sprintf('npp_%d.mat', YEARS(y)));
        d = load(f, 'nppG1dg');
        npp_mo(:,:,:,y) = d.nppG1dg;
    end
    fprintf('  NPP loaded in %.0fs\n', toc(t0));

    fprintf('  Loading %d years of MLD ...\n', NYR); t0 = tic;
    mld_mo = nan(NLAT,NLON,NMO,NYR,'single');
    for y = 1:NYR
        f = fullfile(mld_dir.(pname), sprintf('mld_%d.mat', YEARS(y)));
        d = load(f, 'mldG1deg');
        mld_mo(:,:,:,y) = single(d.mldG1deg);
    end
    fprintf('  MLD loaded in %.0fs\n', toc(t0));

    for s = 1:4
        fprintf('  Season %d/4 (%s) ...\n', s, season_names{s});

        % ── Hemisphere-aware seasonal mean per year ──────────────────────────
        sm_npp = nan(NLAT,NLON,NYR,'single');
        sm_mld = nan(NLAT,NLON,NYR,'single');
        sm_npp(idxN,:,:) = mean(npp_mo(idxN,:,mmN(s,:),:), 3);
        sm_npp(idxS,:,:) = mean(npp_mo(idxS,:,mmS(s,:),:), 3);
        sm_mld(idxN,:,:) = mean(mld_mo(idxN,:,mmN(s,:),:), 3);
        sm_mld(idxS,:,:) = mean(mld_mo(idxS,:,mmS(s,:),:), 3);

        [tp,ts,ms,nv] = cell_ols_trend_pct(sm_npp, t_vec, Sxx, dof, NYR);
        NPP.trend_pct(:,:,s,p)  = tp; NPP.trend_se(:,:,s,p) = ts;
        NPP.mean_state(:,:,s,p) = ms; NPP.n_valid(:,:,s,p)  = nv;

        [tp,ts,ms,nv] = cell_ols_trend_pct(sm_mld, t_vec, Sxx, dof, NYR);
        MLD.trend_pct(:,:,s,p)  = tp; MLD.trend_se(:,:,s,p) = ts;
        MLD.mean_state(:,:,s,p) = ms; MLD.n_valid(:,:,s,p)  = nv;

        fprintf('    valid cells: NPP=%d  MLD=%d  (of %d ocean-grid cells)\n', ...
            sum(NPP.n_valid(:,:,s,p)==NYR,'all'), sum(MLD.n_valid(:,:,s,p)==NYR,'all'), NLAT*NLON);
    end
    clear npp_mo mld_mo
end

save(fullfile(OUT_DIR,'grid_trends_mld_npp.mat'), ...
    'lat','lon','season_names','product_names','NPP','MLD','YEARS','t_vec','-v7.3');
fprintf('\nSaved: %s\n', fullfile(OUT_DIR,'grid_trends_mld_npp.mat'));


%% ── Local function: vectorized per-cell OLS trend, % yr^-1, corrected SE ──
function [trend_pct,trend_se,mean_state,n_valid] = cell_ols_trend_pct(SM, t, Sxx, dof, NYR)
% SM: (nlat,nlon,NYR) seasonal-mean-per-year array
sz = size(SM);
Y = reshape(SM, [sz(1)*sz(2), NYR])';           % NYR x Ncell
valid = all(isfinite(Y),1);                      % 1 x Ncell
n_valid = uint8(sum(isfinite(Y),1));

b = nan(1,size(Y,2)); a = nan(1,size(Y,2));
se_b = nan(1,size(Y,2)); se_a = nan(1,size(Y,2));

Yv = Y(:,valid);
bv = (t' * Yv) / Sxx;                             % 1 x Nvalid
av = mean(Yv,1);
fitted = av + t*bv;                               % NYR x Nvalid
resid  = Yv - fitted;
sigma2 = sum(resid.^2,1) / dof;
se_bv  = sqrt(sigma2 / Sxx);
se_av  = sqrt(sigma2 / NYR);

b(valid) = bv; a(valid) = av; se_b(valid) = se_bv; se_a(valid) = se_av;

trend_pct = b ./ a * 100;
trend_se  = abs(trend_pct) .* sqrt((se_b./b).^2 + (se_a./a).^2);
mean_state = a;

trend_pct  = single(reshape(trend_pct,  sz(1), sz(2)));
trend_se   = single(reshape(trend_se,   sz(1), sz(2)));
mean_state = single(reshape(mean_state, sz(1), sz(2)));
n_valid    = reshape(n_valid, sz(1), sz(2));
end
