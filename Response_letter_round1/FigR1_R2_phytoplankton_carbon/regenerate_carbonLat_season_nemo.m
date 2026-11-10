% Used for: seasonal latitudinal phytoplankton carbon, carbonLat_season_nemo.mat (Fig. R2)
% Original: Manuscript\MLD\manuscript\nature\NCC\R1\code for revision\sensitivty_carbon_mld\regenerate_carbonLat_season_nemo.m (data paths updated to I:\Prj; computation unchanged)
% regenerate_carbonLat_season_nemo.m
%
% Exact replica of Code\npp_cbpm\seasonLat_ot.m's Section 3 (CARBON),
% NEMO only, with corrected I:\ paths (original used stale G:\ paths).
% The pre-existing carbonLat_season_nemo.mat was found to be missing the
% Autumn field entirely (incomplete/stale save) -- this regenerates a
% complete file using the EXACT established methodology: average the 3
% months' RAW GRID values at each (lon,lat) cell first (bathymetry-masked
% via topo>-1000, biome-idx reconstructed), THEN zonally average over
% longitude. This order matters and is NOT equivalent to averaging
% already-zonally-reduced monthly lat-series (a bug found and fixed
% 2026-07-23 in LAT_vs_carbon_mkTrends_allSeasons_nemo.m and
% scatter_npp_vs_carbon_trends_allSeasons_nemo.m, which had used the
% wrong order and produced trends that didn't match the established/
% published values).
%
% Haipeng Zhao
% 2026-07-23

clear; clc;
addpath('D:\Users\zhp19\OneDrive - Indiana State University\Research\MATLAB toolbox\m_map'); % m_topo.mat

OUT_DIR = fileparts(mfilename('fullpath'));
PATH_IN = 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\biomes_carbon_nemo';
PATH_OT = OUT_DIR;
FILENAME_OT = 'carbonLat_season_nemo';

d = load('m_topo.mat');
topo = d.topo';

d = load(fullfile(PATH_IN, '2003', 'carbon_A EQU.mat'));
lon = d.lon;
lat = d.lat;
clear d;

u1 = length(lon);
u2 = length(lat);
carbon = nan(length(lon), length(lat), 20*12);
basin_no = nan(length(lon), length(lat));

datafile = dir([PATH_IN, '\2003\', '*.mat']);
filenum = length(datafile);

basin = 0;
for k = 1:20
    yyyy = num2str(k + 2002);
    for j = 1:filenum
        regfile = datafile(j);
        filename = regfile.name;
        d = load(fullfile(PATH_IN, yyyy, filename));
        for i = 1:length(d.idx)
            [I, J] = ind2sub([u1, u2], d.idx(i));
            if (topo(I, J) > -1000)
                continue
            end
            carbon(I, J, (k-1)*12+1 : k*12) = d.carbon(i, :);
            basin_no(I, J) = basin;
        end
        basin = basin + 1;
    end
    basin = 0;
end

mm = repmat([1:12], 1, 20)';

ind1 = mm==12 | mm==1 | mm==2; % Dec, Jan, Feb
carbon1 = carbon(:,:,ind1);
for i = 1:20
    carbonwi(:,:,i) = nanmean(carbon1(:,:,(i-1)*3+1:(i-1)*3+3), 3);
end

ind2 = mm==3 | mm==4 | mm==5;  % Mar, Apr, May
carbon1 = carbon(:,:,ind2);
for i = 1:20
    carbonsp(:,:,i) = nanmean(carbon1(:,:,(i-1)*3+1:(i-1)*3+3), 3);
end

ind3 = mm==6 | mm==7 | mm==8; % Jun, Jul, Aug
carbon1 = carbon(:,:,ind3);
for i = 1:20
    carbonsu(:,:,i) = nanmean(carbon1(:,:,(i-1)*3+1:(i-1)*3+3), 3);
end

ind4 = mm==9 | mm==10 | mm==11; % Sep, Oct, Nov
carbon1 = carbon(:,:,ind4);
for i = 1:20
    carbonau(:,:,i) = nanmean(carbon1(:,:,(i-1)*3+1:(i-1)*3+3), 3);
end

idxS = lat < 0;

carbonwiS = carbonwi(:,idxS,:);
carbonwi(:,idxS,:) = carbonsu(:,idxS,:);
carbonsu(:,idxS,:) = carbonwiS;

carbonspS = carbonsp(:,idxS,:);
carbonsp(:,idxS,:) = carbonau(:,idxS,:);
carbonau(:,idxS,:) = carbonspS;

a = squeeze(sum(~isnan(carbonwi) .* ~isnan(repmat(basin_no,1,1,20))));
b = squeeze(sum(~isnan(repmat(basin_no,1,1,20))));
c = a ./ b;
carbonwiLat = squeeze(nanmean(carbonwi));
carbonwiLat(c < 0.5) = nan;

a = squeeze(sum(~isnan(carbonsp) .* ~isnan(repmat(basin_no,1,1,20))));
c = a ./ b;
carbonspLat = squeeze(nanmean(carbonsp));
carbonspLat(c < 0.5) = nan;

a = squeeze(sum(~isnan(carbonsu) .* ~isnan(repmat(basin_no,1,1,20))));
c = a ./ b;
carbonsuLat = squeeze(nanmean(carbonsu));
carbonsuLat(c < 0.5) = nan;

a = squeeze(sum(~isnan(carbonau) .* ~isnan(repmat(basin_no,1,1,20))));
c = a ./ b;
carbonauLat = squeeze(nanmean(carbonau));
carbonauLat(c < 0.5) = nan;

save(fullfile(PATH_OT, FILENAME_OT), 'carbonwi','carbonsp','carbonsu','carbonau', ...
    'carbonwiLat','carbonspLat','carbonsuLat','carbonauLat','mm','lat');
fprintf('Saved: %s\n', fullfile(PATH_OT, [FILENAME_OT '.mat']));
