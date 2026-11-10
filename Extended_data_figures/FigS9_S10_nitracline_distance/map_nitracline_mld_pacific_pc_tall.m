% Used for: Fig. S9 (title removed by hand)
% Original: Manuscript\MLD\manuscript\nature\NCC\R1\code for revision\mld_npp_relationship_grid\map_nitracline_mld_pacific_pc_tall.m (data paths updated to I:\Prj; computation unchanged)
% map_nitracline_mld_pacific_pc_tall.m
%
% Same as map_nitracline_mld_pacific_platecarree.m (Pacific-centered,
% -40..40 lat, Plate Carree), but with the y (latitude) axis stretched
% 2x relative to the true 1:1 geographic aspect ratio -- an 80 deg-tall x
% 360 deg-wide band is very flat/hard to read at true scale, so this
% trades strict shape accuracy for on-page legibility (a standard
% compromise for narrow-latitude-band maps). Adjust VSTRETCH below to
% taste (2 = 2x taller than true scale, 1 = true scale/no stretch).
%
% Pacific re-centering note: the underlying zno3/mld grids are already
% stored, after the existing column reorder below, in continuous
% 0.5..359.5 deg-East physical column order (0.5..179.5 = 0-180E,
% 180.5..359.5 = 180E..360E = 180W..0W) -- this was originally labeled
% with a signed -180..180 `lon` vector for the Atlantic-centered map.
% Relabeling the SAME columns with an unsigned, monotonically increasing
% 0.5..359.5 `lon` vector and setting m_proj's longitude range to
% [0.5 359.5] puts x=180 (the physical dateline / mid-Pacific) at the
% center of the map instead of the edges -- no data reordering needed.
%
% Output: fig_mld_nitracline_distance_seasonal_lowmidlat_pacific.png
% (this folder).
%
% Haipeng Zhao
% 2026-09-04

clear; close all; clc;
OUT_DIR = fileparts(mfilename('fullpath'));

% --- Nitracline (nofix) ---
dn = load('I:\Prj\WOA23\nitracline\zno3_woa23_3muM_nofix.mat');
zno3 = dn.zno3(:, [181:360, 1:180], :);   % reorder to continuous 0.5..359.5 E column order
lat = dn.lat(:);
lon = (0.5:1:359.5)';                     % unsigned, Pacific-centered longitude labeling

% --- MLD-NEMO ---
MLD_DIR = 'I:\Prj\MLD_1M_NEMO\one_d03\mld';
YEARS = 2003:2022;
fprintf('Loading %d years of MLD-NEMO ...\n', numel(YEARS));
mld_mo = nan(180,360,12,numel(YEARS),'single');
for y = 1:numel(YEARS)
    d = load(fullfile(MLD_DIR, sprintf('mld_%d.mat', YEARS(y))), 'mldG1deg');
    mld_mo(:,:,:,y) = d.mldG1deg;
end

season_names = {'Spring','Summer','Autumn','Winter'};
mmN = [3 4 5; 6 7 8; 9 10 11; 12 1 2];   % lat >= 0
mmS = [9 10 11; 12 1 2; 3 4 5; 6 7 8];   % lat <  0
idxN = lat >= 0; idxS = lat < 0;

BIOME_DIR = 'I:\Prj\MLD_1M_NEMO\biomes_1M_1dg\mld\2003';
biomes = {'A EQU','EP EQU','WP EQU','NA STPS','SA STPS','NP STPS','SP STPS','IND STPS', ...
    'NA STSS','NP STSS','SO STSS','NA SPSS','NP SPSS','SO SPSS'};
latmask40 = abs(lat) >= 40;   % crop out the +-40-60 bands
VSTRETCH = 2;                 % vertical exaggeration relative to true geographic aspect

fig = figure('Units','centimeters','Position',[2 2 30 12*VSTRETCH*0.6],'Color','w');
t = tiledlayout(fig,2,2,'Padding','compact','TileSpacing','compact');
t.OuterPosition = [0 0.16 1 0.78];  % leave headroom at top (sgtitle) and bottom (colorbar+label)
cmap = crameri('vik', 255);
CLIM = 200;

for s = 1:4
    % Nitracline seasonal mean
    sm_z = nan(180,360,'double');
    sm_z(idxN,:) = mean(zno3(idxN,:,mmN(s,:)), 3, 'omitnan');
    sm_z(idxS,:) = mean(zno3(idxS,:,mmS(s,:)), 3, 'omitnan');

    % MLD seasonal mean (20-year climatology)
    sm_mld_yr = nan(180,360,numel(YEARS),'double');
    sm_mld_yr(idxN,:,:) = mean(mld_mo(idxN,:,mmN(s,:),:), 3, 'omitnan');
    sm_mld_yr(idxS,:,:) = mean(mld_mo(idxS,:,mmS(s,:),:), 3, 'omitnan');
    sm_mld = mean(sm_mld_yr, 3, 'omitnan');

    dist = sm_z - sm_mld;   % positive = MLD shallower than nitracline (light-limited regime)
    distT = dist';          % (360,180) = [lon,lat]

    grid_val = nan(360,180,'single');
    for bi = 1:numel(biomes)
        bidx = load(fullfile(BIOME_DIR, sprintf('mld_%s.mat', biomes{bi})), 'idx');
        idx = round(bidx.idx(:));
        grid_val(idx) = distT(idx);
    end
    grid_val(:, latmask40) = nan;

    ax = nexttile(t);
    m_proj('Equidistant Cylindrical', 'lat', [-40 40], 'lon', [0.5 359.5]);
    hold(ax, 'on');
    m_pcolor(ax, lon, lat, grid_val');
    shading(ax, 'flat');
    m_coast('patch', [0 0 0], 'edgecolor', 'none');
    m_grid('tickdir', 'in', 'linestyle', '--', 'fontsize', 8, 'ytick', -40:20:40, ...
        'xtick', [60 120 180 240 300]);
    colormap(ax, cmap);
    caxis(ax, [-CLIM CLIM]); %#ok<CAXIS>
    set(ax, 'DataAspectRatio', [1 1/VSTRETCH 1]);   % true shape, then stretched VSTRETCH-x taller
    title(ax, season_names{s}, 'FontWeight','normal', 'FontSize', 12);

    if s == 4
        cb = colorbar(ax, 'Location', 'southoutside', 'Position', [0.30 0.09 0.40 0.02]);
        cb.Label.String = 'Nitracline - MLD (m)';
        cb.Label.FontWeight = 'bold';
    end
    fprintf('%s done\n', season_names{s});
end

sgtitle(fig, {'Nitracline - MLD distance (-40 to 40 lat, Pacific-centered, y stretched 2x)', 'brown = MLD shallower than nitracline   |   blue = MLD reaches/exceeds nitracline'}, 'FontWeight','bold', 'FontSize', 12);

out_png = fullfile(OUT_DIR,'fig_mld_nitracline_distance_seasonal_lowmidlat_pacific_tall.png');
exportgraphics(fig, out_png, 'Resolution', 200);
fprintf('Figure saved: %s\n', out_png);
