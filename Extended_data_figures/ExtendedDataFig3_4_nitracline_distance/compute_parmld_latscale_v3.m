% Used for: mixed-layer PAR panel data, parmld_rightpanel_v3.mat (read by the Extended Data Fig. 4 script; that panel is not in the published figure)
% Original: Manuscript\MLD\manuscript\nature\NCC\R1\code for revision\mld_npp_relationship_grid\compute_parmld_latscale_v3.m (data paths updated to I:\Prj; computation unchanged)
% compute_parmld_latscale_v3.m
%
% par_mld computed at the LATITUDINAL scale, not the grid scale.
%
% Difference from v2. v2 evaluated the CbPM light formula per pixel and
% then took the zonal mean:
%       par_mld_lat = < f(PAR, Kd, MLD) >_lon
% v3 takes the zonal mean of each input first and then evaluates the
% formula once per latitude:
%       par_mld_lat =   f( <PAR>_lon, <Kd>_lon, <MLD>_lon )
% These differ because f is nonlinear in MLD (exp(-K*MLD/2)), so this is
% the same Jensen effect that separates the two ways of using a P-E
% curve. Doing it at the latitudinal scale also lets the MLD trend come
% from the paper's own zonal-mean MLD series (mldLat_season.mat, the
% same series the left panel of the figure uses) rather than from
% per-pixel Sen fits.
%
% MLD cases, all from mldLat_season.mat, per latitude per season:
%   climatology  20-year mean of the zonal-mean seasonal MLD
%   2003 / 2022  Theil-Sen (ktaub) endpoints of that same series
% PAR and Kd are held at their climatology throughout, so the spread
% between the 2003 and 2022 curves is the MLD trend acting alone.
%
% Units: the CbPM formula without the daylength division gives
% mol quanta m-2 d-1 directly.
%
% Mask: the paper's convention, all 21 biomes with topo <= -1000 m,
% from topo_biome_mask_nemo.mat.
%
% Output: parmld_rightpanel_v3.mat  (clim, e03, e22, each 180 x 4)

clear; clc;

ROOT    = 'D:\Users\zhp19\OneDrive - Indiana State University\Research\Projects\3D bio model';
OUT_DIR = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(fileparts(OUT_DIR)), 'helper_functions'));                       % austinPetzold_1986


LAMBDA      = [400, 412, 443, 490, 510, 555, 625, 670, 700];
parFraction = [0.0029, 0.0032, 0.0035, 0.0037, 0.0037, 0.0036, 0.0032, 0.0030, 0.0024];

season_names = {'Spring','Summer','Autumn','Winter'};
mmN = [3 4 5; 6 7 8; 9 10 11; 12 1 2];
mmS = [9 10 11; 12 1 2; 3 4 5; 6 7 8];

%% ---------- mask (paper convention) ----------
mk   = load(fullfile(fileparts(OUT_DIR), 'ExtendedDataFig8_grid_cell_scatter', 'topo_biome_mask_nemo.mat'));
mask = logical(mk.valid_lonlat)';        % -> (180 lat, 360 lon)
lat  = mk.lat(:);
idxN = lat >= 0; idxS = lat < 0;

%% ---------- zonal-mean PAR and Kd climatology, per latitude per season ----------
S = load(fullfile(OUT_DIR, 'par_mld_mldonly.mat'), 'PAR_clim', 'KD_clim');
PARm = nan(180, 12); KDm = nan(180, 12);
for m = 1:12
    p = double(S.PAR_clim(:,:,m)); p(~mask) = nan;
    k = double(S.KD_clim(:,:,m));  k(~mask) = nan;
    PARm(:,m) = mean(p, 2, 'omitnan');
    KDm(:,m)  = mean(k, 2, 'omitnan');
end

PAR_lat = nan(180,4); KD_lat = nan(180,4);
for s = 1:4
    PAR_lat(idxN,s) = mean(PARm(idxN, mmN(s,:)), 2, 'omitnan');
    PAR_lat(idxS,s) = mean(PARm(idxS, mmS(s,:)), 2, 'omitnan');
    KD_lat(idxN,s)  = mean(KDm(idxN,  mmN(s,:)), 2, 'omitnan');
    KD_lat(idxS,s)  = mean(KDm(idxS,  mmS(s,:)), 2, 'omitnan');
end

%% ---------- latitudinal MLD: climatology and Sen endpoints ----------
dm = load('I:\Prj\MLD_1M_NEMO\mldLat_season.mat');
mld_vars = {dm.mldspLat, dm.mldsuLat, dm.mldauLat, dm.mldwiLat};

MLD_clim = nan(180,4); MLD_03 = nan(180,4); MLD_22 = nan(180,4);
for s = 1:4
    Y = mld_vars{s};                       % (180 lat, 20 yr)
    MLD_clim(:,s) = mean(Y, 2, 'omitnan');
    for i = 1:180
        y = Y(i,:);
        if sum(isnan(y)) > 2, continue; end
        x = 1:numel(y);
        [~,~,~,~,~,~,~,~,~, pos, ~] = ktaub([x', y'], 0.9999, 1);
        if ~isnan(pos)
            MLD_03(i,s) = pos(1,2);
            MLD_22(i,s) = pos(end,2);
        end
    end
    fprintf('%s latitudinal MLD trend done\n', season_names{s});
end

%% ---------- evaluate the CbPM light formula at the latitudinal scale ----------
clim = parmld_lat(PAR_lat, KD_lat, MLD_clim, LAMBDA, parFraction);
e03  = parmld_lat(PAR_lat, KD_lat, MLD_03,   LAMBDA, parFraction);
e22  = parmld_lat(PAR_lat, KD_lat, MLD_22,   LAMBDA, parFraction);

save(fullfile(OUT_DIR, 'parmld_rightpanel_v3.mat'), 'clim', 'e03', 'e22', 'lat', ...
     'MLD_clim', 'MLD_03', 'MLD_22', 'PAR_lat', 'KD_lat', 'season_names');
fprintf('\nsaved parmld_rightpanel_v3.mat (mol quanta m-2 d-1)\n\n');

%% ---------- 20-40 band summary, cos-lat weighted ----------
band = abs(lat) >= 20 & abs(lat) < 40;
w    = cosd(lat);
fprintf('%-8s %10s %10s %10s %14s %12s\n','season','clim','MLD-2003','MLD-2022','change','MLD 03->22');
for s = 1:4
    v = band & isfinite(clim(:,s)) & isfinite(e03(:,s)) & isfinite(e22(:,s));
    c = sum(clim(v,s).*w(v))/sum(w(v));
    a = sum(e03(v,s).*w(v))/sum(w(v));
    b = sum(e22(v,s).*w(v))/sum(w(v));
    m1 = sum(MLD_03(v,s).*w(v))/sum(w(v));
    m2 = sum(MLD_22(v,s).*w(v))/sum(w(v));
    fprintf('%-8s %10.2f %10.2f %10.2f %7.2f (%+5.1f%%) %5.1f->%5.1f m\n', ...
        season_names{s}, c, a, b, b-a, 100*(b-a)/a, m1, m2);
end

%% ---------- local function ----------
function out = parmld_lat(PAR, KD, MLD, LAMBDA, parFraction)
% Spectrally resolved CbPM mixed-layer light, evaluated on (180 lat x 4
% season) arrays. No daylength division, so the result is the daily
% integral in mol quanta m-2 d-1.
out = nan(size(PAR));
for s = 1:size(PAR,2)
    p = PAR(:,s); k = KD(:,s); d = MLD(:,s);
    acc = zeros(size(p)); Ezprev = [];
    for q = 1:numel(LAMBDA)
        Kq  = austinPetzold_1986(LAMBDA(q), k);
        Ezq = p .* parFraction(q) .* 0.975 .* exp(-Kq .* d / 2);
        if q > 1
            acc = acc + (LAMBDA(q)-LAMBDA(q-1)) .* (Ezq + Ezprev)/2;
        end
        Ezprev = Ezq;
    end
    out(:,s) = acc;
end
end
