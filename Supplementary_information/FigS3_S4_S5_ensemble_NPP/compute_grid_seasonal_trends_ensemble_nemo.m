% Used for: per-grid-cell seasonal trends of ensemble NPP and MLD-NEMO (SI Fig. S4)
% Original: Manuscript\MLD\manuscript\nature\NCC\R1\code for revision\mld_npp_relationship_grid\compute_grid_seasonal_trends_ensemble_nemo.m (data paths updated to I:\Prj; computation unchanged)
% compute_grid_seasonal_trends_ensemble_nemo.m
%
% Same grid-cell (1x1 deg) seasonal trend methodology as
% compute_grid_seasonal_trends.m, but NPP = VGPM+CAFE+CbPM ensemble mean
% instead of CbPM alone, paired only with MLD-NEMO (not ECCO/HYCOM).
%
% Inputs (read-only):
%   NPP (VGPM+CAFE+CbPM ensemble): I:\Prj\npp_ensemble_1M_1dg\grid_npp_nemo_1M_1dg_mat\npp_YYYY.mat  (nppG1dg)
%   MLD (NEMO)                   : I:\Prj\MLD_1M_NEMO\one_d03\mld\mld_YYYY.mat                         (mldG1deg)
% CORRECTED 2026-07-25: was I:\Prj\MLD_1M_NEMO\one\ (stale, superseded by
% one_d03\mld\) -- same fix/rationale as compute_grid_seasonal_trends.m.
%
% Output: grid_trends_mld_npp_ensemble_nemo.mat (this folder) -- lat, lon,
% season_names, NPP/MLD trend_pct / trend_se / mean_state / n_valid grids,
% size (180,360,4) -- single product (NEMO), no product dimension.
%
% Haipeng Zhao
% 2026-07-22

clear; close all; clc;
OUT_DIR = fullfile(fileparts(mfilename('fullpath')));

YEARS   = 2003:2022;
NYR     = numel(YEARS);
NLAT    = 180; NLON = 360; NMO = 12;

npp_dir = 'I:\Prj\npp_ensemble_1M_1dg\grid_npp_nemo_1M_1dg_mat';
mld_dir = 'I:\Prj\MLD_1M_NEMO\one_d03\mld';  % CORRECTED 2026-07-25, was 'I:\Prj\MLD_1M_NEMO\one'

season_names = {'Spring','Summer','Autumn','Winter'};
mmN = [3 4 5; 6 7 8; 9 10 11; 12 1 2];   % lat >= 0
mmS = [9 10 11; 12 1 2; 3 4 5; 6 7 8];   % lat <  0

% Reference grid
ref = load(fullfile(npp_dir, sprintf('npp_%d.mat', YEARS(1))), 'lat','lon');
lat = double(ref.lat(:));   % (180,1), -89.5 .. 89.5
lon = double(ref.lon(:));   % (360,1),   0.5 .. 359.5
assert(numel(lat)==NLAT && numel(lon)==NLON);
idxN = lat >= 0;
idxS = lat <  0;

t_vec = (-9.5:9.5)';                 % centered year index, NYR=20
Sxx   = sum(t_vec.^2);
dof   = NYR - 2;

NPP.trend_pct  = nan(NLAT,NLON,4,'single');
NPP.trend_se   = nan(NLAT,NLON,4,'single');
NPP.mean_state = nan(NLAT,NLON,4,'single');
NPP.n_valid    = zeros(NLAT,NLON,4,'uint8');
MLD = NPP;

fprintf('Loading %d years of ensemble NPP (NEMO-paired) ...\n', NYR); t0 = tic;
npp_mo = nan(NLAT,NLON,NMO,NYR,'single');
for y = 1:NYR
    f = fullfile(npp_dir, sprintf('npp_%d.mat', YEARS(y)));
    d = load(f, 'nppG1dg');
    npp_mo(:,:,:,y) = d.nppG1dg;
end
fprintf('NPP loaded in %.0fs\n', toc(t0));

fprintf('Loading %d years of MLD-NEMO ...\n', NYR); t0 = tic;
mld_mo = nan(NLAT,NLON,NMO,NYR,'single');
for y = 1:NYR
    f = fullfile(mld_dir, sprintf('mld_%d.mat', YEARS(y)));
    d = load(f, 'mldG1deg');
    mld_mo(:,:,:,y) = single(d.mldG1deg);
end
fprintf('MLD loaded in %.0fs\n', toc(t0));

for s = 1:4
    fprintf('Season %d/4 (%s) ...\n', s, season_names{s});

    sm_npp = nan(NLAT,NLON,NYR,'single');
    sm_mld = nan(NLAT,NLON,NYR,'single');
    sm_npp(idxN,:,:) = mean(npp_mo(idxN,:,mmN(s,:),:), 3);
    sm_npp(idxS,:,:) = mean(npp_mo(idxS,:,mmS(s,:),:), 3);
    sm_mld(idxN,:,:) = mean(mld_mo(idxN,:,mmN(s,:),:), 3);
    sm_mld(idxS,:,:) = mean(mld_mo(idxS,:,mmS(s,:),:), 3);

    [tp,ts,ms,nv] = cell_ols_trend_pct(sm_npp, t_vec, Sxx, dof, NYR);
    NPP.trend_pct(:,:,s)  = tp; NPP.trend_se(:,:,s) = ts;
    NPP.mean_state(:,:,s) = ms; NPP.n_valid(:,:,s)  = nv;

    [tp,ts,ms,nv] = cell_ols_trend_pct(sm_mld, t_vec, Sxx, dof, NYR);
    MLD.trend_pct(:,:,s)  = tp; MLD.trend_se(:,:,s) = ts;
    MLD.mean_state(:,:,s) = ms; MLD.n_valid(:,:,s)  = nv;

    fprintf('  valid cells: NPP=%d  MLD=%d  (of %d ocean-grid cells)\n', ...
        sum(NPP.n_valid(:,:,s)==NYR,'all'), sum(MLD.n_valid(:,:,s)==NYR,'all'), NLAT*NLON);
end

save(fullfile(OUT_DIR,'grid_trends_mld_npp_ensemble_nemo.mat'), ...
    'lat','lon','season_names','NPP','MLD','YEARS','t_vec','-v7.3');
fprintf('\nSaved: %s\n', fullfile(OUT_DIR,'grid_trends_mld_npp_ensemble_nemo.mat'));


%% -- Local function: vectorized per-cell OLS trend, % yr^-1, corrected SE --
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
