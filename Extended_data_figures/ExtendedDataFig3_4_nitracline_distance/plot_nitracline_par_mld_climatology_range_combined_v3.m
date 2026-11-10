% Used for: Extended Data Fig. 4 (left panel; the PAR panel was removed by hand)
% Original: Manuscript\MLD\manuscript\nature\NCC\R1\code for revision\mld_npp_relationship_grid\plot_nitracline_par_mld_climatology_range_combined_v3.m (data paths updated to I:\Prj; computation unchanged)
% plot_nitracline_par_mld_climatology_range_combined.m
%
% Two-panel, -40..40 lat, merged-seasons (IPCC AR4 SPM-style) figure:
%   left  = zonal-mean Nitracline - MLD distance
%   right = zonal-mean par_mld (CbPM's mixed-layer-averaged light climate)
%
% Both panels use the SAME convention:
%   solid line  = 2003-2022 climatology (plain 20-yr mean per latitude row)
%   shaded band = range spanned by the paper's own per-latitude ktaub
%                 (Theil-Sen) 2003/2022 MLD-trend endpoints -- i.e. how
%                 much that quantity moved due to the MLD trend alone
%
% Left panel: nitracline held at its WOA23 climatology (unaffected by the
% MLD-trend question), MLD from mldLat_season.mat (zonal-mean-first,
% paper's own method) -- same construction as
% plot_nitracline_mld_latavg_lowmidlat_merged.m.
%
% MASK HARMONISED 2026-09-12: the nitracline zonal mean used a 14-biome
% list with no topo filter, while the MLD series (mldLat_season.mat) and
% the whole right panel use the paper's convention (all 21 biomes, topo
% <= -1000 m). The two terms of "nitracline - MLD" were therefore averaged
% over different cell sets. Now both use topo_biome_mask_nemo.mat. Effect
% is small (the old mask was a strict superset within 40S-40N, 428 extra
% shelf cells of ~18000, shifting zonal-mean nitracline depth by 0.3-0.6 m)
% and the shaded band is unchanged exactly, since the nitracline term
% cancels in dist2022 - dist2003.
%
% Right panel: par_mld already reflects each year's actual MLD (it's a
% direct per-pixel model OUTPUT, not a "hold X fixed" construction), from
% par_mldLat_season.mat -- climatology = 20-yr mean, band = 2003/2022
% ktaub trend endpoints, same as the left panel's convention.
%
% Single shared legend (4 seasons + shading swatch) at the bottom.
%
% 2026-09-11 revisions (styling only, no change to any computed value):
%   - overall figure title removed (explained in the caption instead)
%   - latitude ticks labelled 40S..0..40N on both panels
%   - larger fonts throughout
%   - shading legend entry reworded to "Change due to MLD trend (2003 to 2022)"
%   - left-panel arrays cached to nitracline_distance_leftpanel_v3.mat
%     after each full computation, so later restyling runs work without
%     the portable drive (I:) connected
%   - ktaub path resolved on either machine (desktop D:\Users\zhp19\...,
%     laptop E:\...)
%
% Output: fig_nitracline_par_mld_climatology_range_combined_v3.png (this folder).
%
% Haipeng Zhao
% 2026-09-04, restyled 2026-09-11

clear; close all; clc;
OUT_DIR = fileparts(mfilename('fullpath'));

season_names = {'Spring','Summer','Autumn','Winter'};
% Okabe-Ito-derived, colorblind-safe categorical palette (same validated
% set used elsewhere in this project's response figures).
season_colors = [0.784 0.498 0.000; 0.000 0.620 0.451; 0.30 0.30 0.30; 0.000 0.447 0.698];
mmN = [3 4 5; 6 7 8; 9 10 11; 12 1 2];
mmS = [9 10 11; 12 1 2; 3 4 5; 6 7 8];

%% ---------- LEFT: Nitracline - MLD ----------
NITRACLINE_FILE = 'I:\Prj\WOA23\nitracline\zno3_woa23_3muM_nofix.mat';
LEFT_CACHE = fullfile(OUT_DIR, 'nitracline_distance_leftpanel_v3.mat');

if exist(NITRACLINE_FILE, 'file')
    for kp = {'D:\Users\zhp19\OneDrive - Indiana State University\Research\MATLAB toolbox\ktaub', ...
              'E:\OneDrive - Indiana State University\Research\MATLAB toolbox\ktaub'}
        if exist(kp{1}, 'dir'), addpath(kp{1}); end
    end

    dn = load(NITRACLINE_FILE);
    zno3 = dn.zno3(:, [181:360, 1:180], :);
    lat = dn.lat(:);
    idxN = lat >= 0; idxS = lat < 0;

    mk = load(fullfile(OUT_DIR, '..', 'ExtendedDataFig8_grid_cell_scatter', 'topo_biome_mask_nemo.mat'), 'valid_lonlat');
    mask = logical(mk.valid_lonlat);   % (360 lon, 180 lat)

    nitracline_lat = nan(180, 4);
    for s = 1:4
        sm_z = nan(180,360,'double');
        sm_z(idxN,:) = mean(zno3(idxN,:,mmN(s,:)), 3, 'omitnan');
        sm_z(idxS,:) = mean(zno3(idxS,:,mmS(s,:)), 3, 'omitnan');
        sm_zT = sm_z';
        sm_zT(~mask) = nan;
        nitracline_lat(:, s) = mean(sm_zT, 1, 'omitnan')';
    end

    dm = load('I:\Prj\MLD_1M_NEMO\mldLat_season.mat');
    mld_lat_vars = {dm.mldspLat, dm.mldsuLat, dm.mldauLat, dm.mldwiLat};

    mld2003_lat = nan(180, 4);
    mld2022_lat = nan(180, 4);
    for s = 1:4
        Y = mld_lat_vars{s};
        for i = 1:180
            y = Y(i, :);
            if sum(isnan(y)) > 2, continue; end
            x = 1:length(y);
            [~,~,~,~,~,~,~,~,~, plotofslope, ~] = ktaub([x', y'], 0.9999, 1);
            if ~isnan(plotofslope)
                mld2003_lat(i, s) = plotofslope(1, 2);
                mld2022_lat(i, s) = plotofslope(end, 2);
            end
        end
        fprintf('%s MLD lat-trend done\n', season_names{s});
    end

    mld_clim_lat = nan(180, 4);
    for s = 1:4
        mld_clim_lat(:, s) = mean(mld_lat_vars{s}, 2, 'omitnan');
    end

    dist2003_lat = nitracline_lat - mld2003_lat;
    dist2022_lat = nitracline_lat - mld2022_lat;
    dist_clim_lat = nitracline_lat - mld_clim_lat;

    save(LEFT_CACHE, 'lat', 'dist2003_lat', 'dist2022_lat', 'dist_clim_lat');
    fprintf('Left-panel arrays cached to %s\n', LEFT_CACHE);
elseif exist(LEFT_CACHE, 'file')
    load(LEFT_CACHE, 'lat', 'dist2003_lat', 'dist2022_lat', 'dist_clim_lat');
    fprintf('Portable drive not found; left panel loaded from cache %s\n', LEFT_CACHE);
else
    error(['Left-panel inputs unavailable: %s not found and no cache at %s. ' ...
           'Connect the portable drive (I:) and run once to build the cache.'], ...
           NITRACLINE_FILE, LEFT_CACHE);
end

%% ---------- RIGHT: average MLD PAR  (v3, latitudinal scale) ----------
% Replaces the old par_mldLat_season.mat computation. Three corrections,
% all applied in export_parmld_rightpanel.py and loaded here:
%   MASK  the paper's convention (all 21 biomes, topo <= -1000 m) from
%         topo_biome_mask_nemo.mat, verified against
%         Code/npp_cbpm/monthlyLat_nemo_ot.m. The old file used a
%         14-biome mask with no topo filter.
%   UNITS par_mld is stored in mol quanta m-2 h-1. It is multiplied by
%         daylength (hours) per month and per latitude BEFORE any
%         averaging, giving mol quanta m-2 d-1. The v1 figure labelled
%         the axis d-1 while plotting h-1 values.
%   RANGE the band is the MLD-ONLY counterfactual (PAR and Kd pinned at
%         climatology, only MLD on its Sen line).
%   SCALE v3 evaluates the CbPM light formula AT THE LATITUDINAL SCALE:
%         PAR, Kd and MLD are zonally averaged first, then the formula is
%         applied once per latitude, and the MLD trend comes from the
%         paper's own zonal-mean MLD series (mldLat_season.mat, the same
%         series the left panel uses). v2 evaluated per pixel and then
%         averaged. See compute_parmld_latscale_v3.m.
dp = load(fullfile(OUT_DIR, 'parmld_rightpanel_v3.mat'));
lat_p        = dp.lat(:);
par_clim_lat = dp.clim;    % (180 x 4) mol quanta m-2 d-1
par2003_lat  = dp.e03;
par2022_lat  = dp.e22;
disp('right panel loaded from parmld_rightpanel_v3.mat (mol quanta m-2 d-1)');

%% ---------- Plot: 1x2, seasons merged, climatology line + trend-range band ----------
FS_AX = 14;   % tick labels and axis labels
FS_TI = 15;   % panel titles
FS_LG = 14;   % legend
lat_ticks = -40:10:40;
lat_tick_labels = arrayfun(@lat_label, lat_ticks, 'UniformOutput', false);

latwin = abs(lat) <= 40;
lat_w = lat(latwin);
latwin_p = abs(lat_p) <= 40;
lat_wp = lat_p(latwin_p);

fig = figure('Units','centimeters','Position',[2 2 28 17],'Color','w');
t = tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'loose');

% --- left panel ---
ax1 = nexttile(t);
hold(ax1, 'on');
h1 = gobjects(4,1);
for s = 1:4
    x1 = dist2003_lat(latwin, s);
    x2 = dist2022_lat(latwin, s);
    xc = dist_clim_lat(latwin, s);
    valid = ~isnan(x1) & ~isnan(x2) & ~isnan(xc);
    fill(ax1, [x1(valid); flipud(x2(valid))], [lat_w(valid); flipud(lat_w(valid))], ...
        season_colors(s,:), 'FaceAlpha', 0.55, 'EdgeColor', season_colors(s,:), ...
        'EdgeAlpha', 0.55, 'LineWidth', 0.5, 'HandleVisibility', 'off');
    h1(s) = plot(ax1, xc(valid), lat_w(valid), '-', 'Color', season_colors(s,:), 'LineWidth', 1.6);
end
ylim(ax1, [-40 40]);
xlim(ax1, [0 200]);   % only nitracline deeper than the MLD
ylabel(ax1, 'Latitude');
xlabel(ax1, 'Nitracline - MLD (m)');
title(ax1, 'Nitracline - MLD distance', 'FontWeight', 'normal', 'FontSize', FS_TI);
style_lat_axis(ax1, lat_ticks, lat_tick_labels, FS_AX);

% --- right panel ---
ax2 = nexttile(t);
hold(ax2, 'on');
for s = 1:4
    x1 = par2003_lat(latwin_p, s);
    x2 = par2022_lat(latwin_p, s);
    xc = par_clim_lat(latwin_p, s);
    valid = ~isnan(x1) & ~isnan(x2) & ~isnan(xc);
    fill(ax2, [x1(valid); flipud(x2(valid))], [lat_wp(valid); flipud(lat_wp(valid))], ...
        season_colors(s,:), 'FaceAlpha', 0.55, 'EdgeColor', season_colors(s,:), ...
        'EdgeAlpha', 0.55, 'LineWidth', 0.5, 'HandleVisibility', 'off');
    plot(ax2, xc(valid), lat_wp(valid), '-', 'Color', season_colors(s,:), 'LineWidth', 1.6);
end
ylim(ax2, [-40 40]);
xlim(ax2, [0 30]);
xlabel(ax2, 'Mixed-layer average PAR (mol quanta m^{-2} d^{-1})');
title(ax2, 'Mixed-layer average PAR', 'FontWeight', 'normal', 'FontSize', FS_TI);
style_lat_axis(ax2, lat_ticks, lat_tick_labels, FS_AX);

% one generic legend swatch explaining the shading, shared by both panels
hband = fill(ax1, nan(4,1), nan(4,1), [0.5 0.5 0.5], 'FaceAlpha', 0.55, 'EdgeColor', 'none');

lg = legend(ax1, [h1; hband], [season_names, {'Change due to MLD trend (2003 to 2022)'}], ...
    'Orientation', 'horizontal', 'Box', 'off', 'NumColumns', 5, 'FontSize', FS_LG);
lg.Layout.Tile = 'south';

out_png = fullfile(OUT_DIR, 'fig_nitracline_par_mld_climatology_range_combined_v3.png');
exportgraphics(fig, out_png, 'Resolution', 200);
fprintf('Figure saved: %s\n', out_png);

function style_lat_axis(ax, ticks, labels, fs)
    yticks(ax, ticks);
    yticklabels(ax, labels);
    ax.FontSize = fs;
    ax.LabelFontSizeMultiplier = 1.0;
    box(ax, 'on');
    ax.Toolbar.Visible = 'off';
end

function s = lat_label(v)
    if v < 0
        s = sprintf('%d%cS', -v, char(176));
    elseif v > 0
        s = sprintf('%d%cN', v, char(176));
    else
        s = sprintf('0%c', char(176));
    end
end
