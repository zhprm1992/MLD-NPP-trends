% Used for: nitracline-MLD distance per grid cell and year, plot_distance_r9style_range.mat (Fig. R3)
% Original: Manuscript\MLD\manuscript\nature\NCC\R2\code for revision\distance_figure\plot_distance_r9style_range.m (data paths updated to I:\Prj; computation unchanged)
% plot_distance_r9style_range.m  (version with the full year-to-year range, and an anomaly version; HZ 2026-11-06)
% based on plot_distance_r9style.m
% Nitracline - MLD distance, 40S-40N, in the style of R1 Extended Data R9
% (R1\code for revision\mld_npp_relationship_grid\plot_nitracline_par_mld_climatology_range_combined_v3.m, left panel).
% Pre-specified in revision_log_R2.md (2026-11-06). WOA23 3 umol/kg nofix nitracline, MLD-NEMO, R9 mask.
% Per grid cell x season (seasonal mean MLD per year, seasonal mean climatological nitracline), then zonal mean:
%   solid line  = climatological distance
%   dark band   = distance on the MLD trend line in 2003 and 2022 (change due to the MLD trend), Theil-Sen
%   light band  = climatology +- 1 SD of the detrended year-to-year distance
%   dashed line = climatology minus the deepest-year excursion (distance in the deepest year, trend removed)
% READ ONLY on inputs. Haipeng Zhao / Claude, 2026-11-06 (R2)

clear; close all; clc;
OUT = fileparts(mfilename('fullpath'));
R1G = fullfile(fileparts(fileparts(OUT)), 'Extended_data_figures', 'ExtendedDataFig8_grid_cell_scatter');
diary(fullfile(OUT, 'plot_distance_r9style_range_log.txt')); diary on
seasName = {'Spring','Summer','Autumn','Winter'};
season_colors = [0.784 0.498 0.000; 0.000 0.620 0.451; 0.30 0.30 0.30; 0.000 0.447 0.698];   % as R9
mmN = [3 4 5; 6 7 8; 9 10 11; 12 1 2]; mmS = [9 10 11; 12 1 2; 3 4 5; 6 7 8];
dn = load('I:\Prj\WOA23\nitracline\zno3_woa23_3muM_nofix.mat'); zno3 = double(dn.zno3(:, [181:360, 1:180], :)); lat = double(dn.lat(:));
mk = load(fullfile(R1G, 'topo_biome_mask_nemo.mat'), 'valid_lonlat'); mask = logical(mk.valid_lonlat)';   % (lat, lon)
NY = 20; M = nan(180, 360, 12, NY);
for y = 1:NY, d = load(['I:\Prj\MLD_1M_NEMO\one_d03\mld\mld_' num2str(2002 + y) '.mat'], 'mldG1deg'); M(:,:,:,y) = d.mldG1deg; end
isN = lat >= 0; t = 1:NY; [i1, j1] = find(triu(true(NY), 1));
Zc = nan(180, 6, 4);                               % per lat row: clim, d2003, d2022, SD, deepest, shallowest
for s = 1:4
    Ms = nan(180, 360, NY); zs = nan(180, 360);
    Ms(isN,:,:) = squeeze(mean(M(isN,:,mmN(s,:),:), 3)); Ms(~isN,:,:) = squeeze(mean(M(~isN,:,mmS(s,:),:), 3));
    zs(isN,:) = mean(zno3(isN,:,mmN(s,:)), 3, 'omitnan'); zs(~isN,:) = mean(zno3(~isN,:,mmS(s,:)), 3, 'omitnan');
    D = reshape(zs, [], 1) - reshape(Ms, [], NY);  % distance per cell and year
    v = all(isfinite(D), 2) & reshape(mask, [], 1) & abs(repmat(lat, 360, 1)) <= 40;
    sl = nan(size(D,1), 1); b0 = sl;
    sl(v) = median((D(v, j1) - D(v, i1)) ./ (t(j1) - t(i1)), 2); b0(v) = median(D(v,:) - sl(v) .* t, 2);
    res = D - (b0 + sl .* t);
    Q = nan(size(D,1), 6);
    Q(v,:) = [mean(D(v,:), 2), b0(v) + sl(v), b0(v) + sl(v) * NY, std(res(v,:), 0, 2), mean(D(v,:), 2) + min(res(v,:), [], 2), mean(D(v,:), 2) + max(res(v,:), [], 2)];
    for j = 1:6, Zc(:, j, s) = mean(reshape(Q(:,j), 180, 360), 2, 'omitnan'); end
    w = abs(lat) < 40 & abs(lat) > 10;
    fprintf('%s, 10-40 both hemispheres (zonal means): |trend change| %.1f m, SD %.1f m, deepest-year excursion %.1f m\n', seasName{s}, ...
        median(abs(Zc(w,3,s) - Zc(w,2,s)), 'omitnan'), median(Zc(w,4,s), 'omitnan'), median(Zc(w,1,s) - Zc(w,5,s), 'omitnan'));
end
save(fullfile(OUT, 'plot_distance_r9style_range.mat'), 'Zc', 'lat');

FS_AX = 12; FS_TI = 13;
lat_ticks = -40:10:40; lat_tick_labels = arrayfun(@lat_label, lat_ticks, 'UniformOutput', false);
w = abs(lat) <= 40; L = lat(w); g = [0.45 0.45 0.45];
for ver = 1:2
    fig = figure('Units', 'centimeters', 'Position', [2 2 32 15], 'Color', 'w');
    tl = tiledlayout(fig, 1, 4, 'Padding', 'compact', 'TileSpacing', 'compact');
    for s = 1:4
        ax = nexttile(tl); hold(ax, 'on'); c = season_colors(s,:);
        q = Zc(w,:,s); ok = all(isfinite(q), 2); q = q(ok,:); Ly = L(ok);
        if ver == 2, q = q - q(:,1); end                                   % anomaly relative to climatology
        fill(ax, [q(:,5); flipud(q(:,6))], [Ly; flipud(Ly)], c, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        fill(ax, [q(:,2); flipud(q(:,3))], [Ly; flipud(Ly)], c, 'FaceAlpha', 0.8, 'EdgeColor', c, 'EdgeAlpha', 0.8, 'LineWidth', 0.5);
        if ver == 1, plot(ax, q(:,1), Ly, '-', 'Color', c, 'LineWidth', 1.6); else, xline(ax, 0, 'k-', 'LineWidth', 0.6); end
        plot(ax, q(:,5), Ly, '--', 'Color', c, 'LineWidth', 1.3);
        ylim(ax, [-40 40]); yticks(ax, lat_ticks); box(ax, 'on'); ax.FontSize = FS_AX; ax.Toolbar.Visible = 'off';
        if ver == 1, xlim(ax, [0 200]); xlabel(ax, 'Nitracline - MLD (m)'); else, xlim(ax, [-50 50]); xlabel(ax, {'Change in nitracline - MLD (m)', '(negative = closer to the nitracline)'}); end
        title(ax, seasName{s}, 'FontWeight', 'normal', 'FontSize', FS_TI);
        if s == 1, yticklabels(ax, lat_tick_labels); ylabel(ax, 'Latitude'); else, yticklabels(ax, {}); end
        text(ax, 0, 1.005, ['(' char('a' + s - 1) ')'], 'Units', 'normalized', 'VerticalAlignment', 'bottom', 'FontWeight', 'bold', 'FontSize', 12)
    end
    l2 = fill(nan(3,1), nan(3,1), g, 'FaceAlpha', 0.8, 'EdgeColor', 'none'); l3 = fill(nan(3,1), nan(3,1), g, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    l4 = plot(nan, nan, '--', 'Color', g, 'LineWidth', 1.3);
    if ver == 1
        l1 = plot(nan, nan, '-', 'Color', g, 'LineWidth', 1.6);
        lg = legend([l1 l2 l3 l4], {'Climatology (2003 to 2022)', 'Change due to MLD trend (2003 to 2022)', 'Year-to-year range', 'Deepest year'}, 'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 12);
        fn = 'fig_distance_r9style_range.png';
    else
        lg = legend([l2 l3 l4], {'Change due to MLD trend (2003 to 2022)', 'Year-to-year range', 'Deepest year'}, 'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 12);
        fn = 'fig_distance_anomaly_range.png';
    end
    lg.Layout.Tile = 'south';
    exportgraphics(fig, fullfile(OUT, fn), 'Resolution', 300);
end
diary off

function s = lat_label(v)
    if v < 0, s = sprintf('%d%cS', -v, char(176)); elseif v > 0, s = sprintf('%d%cN', v, char(176)); else, s = sprintf('0%c', char(176)); end
end
