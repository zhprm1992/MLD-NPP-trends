% Used for: seasonal latitudinal N2, bvfLat_season_weighted.mat (SI Table S1)
% Original: Manuscript\MLD\manuscript\nature\NCC\R1\code for revision\correlation_bvf_npp\latSeason_bvf200_weighted_ot.m (data paths updated to I:\Prj; computation unchanged)
% latSeason_bvf200_weighted_ot.m
%
% Stage 3 of 3: seasonal latitudinal average of the depth-weighted 0-200m
% BVF, mirroring Code\mld_NEMO\latSeason_ot.m EXACTLY (topo>-1000m depth
% mask, basin_no bookkeeping, DJF/MAM/JJA/SON seasonal averaging with
% hemisphere season-swap, 50% per-lat-band completeness threshold) -- NOT
% the simpler latSeason_N2_ot.m (which has neither the topo mask nor the
% coverage threshold). Per user request 2026-07-29: match the same
% methodology as the established MLD/NPP pipeline for full consistency,
% rather than replicate the older, simpler N2-specific script.
%
% Reads from biomes_1M_1dg\bvf200_weighted\<year>\bvf200_<biome>.mat
% (fields bvf200, lat, lon, idx), built by biomes_1M_bvf200_weighted_ot.m.
%
% Output: bvfLat_season_weighted.mat (this folder) -- bvf{wi,sp,su,au},
% bvf{wi,sp,su,au}Lat, mm, lat (same field-naming convention as the
% original bvfLat_season.mat, so it's a drop-in replacement for
% correlation_bvf_npp_nemo_run.m).
%
% Haipeng Zhao
% 2026-07-29

clear; clc;
addpath('D:\Users\zhp19\OneDrive - Indiana State University\Research\MATLAB toolbox\m_map'); % m_topo.mat
OUT_DIR = fileparts(mfilename('fullpath'));

PATH_IN = 'I:\Prj\N2_1M_NEMO\biomes_1M_1dg\bvf200_weighted';

d = load('m_topo.mat');
topo = d.topo';

ref = load(fullfile(PATH_IN, '2003', 'bvf200_A EQU.mat'), 'lat','lon');
lon = ref.lon; lat = ref.lat;

u1 = length(lon); u2 = length(lat);
bvf = nan(u1, u2, 20*12);
basin_no = nan(u1, u2);

datafile = dir(fullfile(PATH_IN, '2003', '*.mat'));
filenum = length(datafile);
fprintf('Found %d biome files\n', filenum);

basin = 0;
for k = 1:20
    yyyy = num2str(k + 2002);
    for j = 1:filenum
        regfile = datafile(j);
        filename = regfile.name;
        disp(filename)

        d = load(fullfile(PATH_IN, yyyy, filename));
        for i = 1:length(d.idx)
            [I, J] = ind2sub([u1, u2], d.idx(i));
            if (topo(I, J) > -1000)
                continue
            end
            bvf(I, J, (k-1)*12+1 : k*12) = d.bvf200(i, :);
            basin_no(I, J) = basin;
        end
        basin = basin + 1;
    end
    basin = 0;
end

mm = repmat(1:12, 1, 20)';

idxN = lat > 0;
idxS = lat < 0;

ind1 = mm == 12 | mm == 1 | mm == 2; % Dec, Jan, Feb
bvf1 = bvf(:, :, ind1);
for i = 1:20
    bvfwi(:, :, i) = nanmean(bvf1(:, :, (i-1)*3+1:(i-1)*3+3), 3);
end

ind2 = mm == 3 | mm == 4 | mm == 5;  % Mar, April, May
bvf1 = bvf(:, :, ind2);
for i = 1:20
    bvfsp(:, :, i) = nanmean(bvf1(:, :, (i-1)*3+1:(i-1)*3+3), 3);
end

ind3 = mm == 6 | mm == 7 | mm == 8; % June, July, Aug
bvf1 = bvf(:, :, ind3);
for i = 1:20
    bvfsu(:, :, i) = nanmean(bvf1(:, :, (i-1)*3+1:(i-1)*3+3), 3);
end

ind4 = mm == 9 | mm == 10 | mm == 11; % Sep, Oct, Nov
bvf1 = bvf(:, :, ind4);
for i = 1:20
    bvfau(:, :, i) = nanmean(bvf1(:, :, (i-1)*3+1:(i-1)*3+3), 3);
end

bvfwiS = bvfwi(:, idxS, :);
bvfwi(:, idxS, :) = bvfsu(:, idxS, :);
bvfsu(:, idxS, :) = bvfwiS;

bvfspS = bvfsp(:, idxS, :);
bvfsp(:, idxS, :) = bvfau(:, idxS, :);
bvfau(:, idxS, :) = bvfspS;

a = squeeze(sum(~isnan(bvfwi) .* ~isnan(repmat(basin_no, 1, 1, 20))));
b = squeeze(sum(~isnan(repmat(basin_no, 1, 1, 20))));
c = a ./ b;
bvfwiLat = squeeze(nanmean(bvfwi));
bvfwiLat(c < 0.5) = nan;

a = squeeze(sum(~isnan(bvfsp) .* ~isnan(repmat(basin_no, 1, 1, 20))));
c = a ./ b;
bvfspLat = squeeze(nanmean(bvfsp));
bvfspLat(c < 0.5) = nan;

a = squeeze(sum(~isnan(bvfsu) .* ~isnan(repmat(basin_no, 1, 1, 20))));
c = a ./ b;
bvfsuLat = squeeze(nanmean(bvfsu));
bvfsuLat(c < 0.5) = nan;

a = squeeze(sum(~isnan(bvfau) .* ~isnan(repmat(basin_no, 1, 1, 20))));
c = a ./ b;
bvfauLat = squeeze(nanmean(bvfau));
bvfauLat(c < 0.5) = nan;

save(fullfile(OUT_DIR, 'bvfLat_season_weighted'), 'bvfwi', 'bvfsp', 'bvfsu', 'bvfau', 'bvfwiLat', 'bvfspLat', 'bvfsuLat', 'bvfauLat', 'mm', 'lat')
fprintf('Saved: %s\n', fullfile(OUT_DIR, 'bvfLat_season_weighted.mat'));
