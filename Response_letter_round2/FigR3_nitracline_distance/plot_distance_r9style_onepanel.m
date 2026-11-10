% Used for: round-2 letter Fig. R3
% Original: Manuscript\MLD\manuscript\nature\NCC\R2\code for revision\distance_figure\plot_distance_r9style_onepanel.m (data paths updated to I:\Prj; computation unchanged)
% plot_distance_r9style_onepanel.m  Figure 10c: nitracline - MLD distance, 40S-40N, all four seasons in one panel.
% Line = climatology (Zc from plot_distance_r9style_range.mat); dark band = distance on the MLD trend line in 2003 and 2022
% (change due to the trend). Light band = range of the ACTUAL yearly distance (trend not removed), as averages:
% left edge = average distance in the years closer to the nitracline than the climatology, right edge = average distance in the
% years farther from it. Per grid cell and season (seasonal mean MLD per year, climatological nitracline), same cells as Zc
% (R9 mask, finite in all years, |lat| <= 40), then zonal mean. HZ decisions 2026-11-07: averages not max-min; actual distance.
% READ ONLY on inputs. Haipeng Zhao / Claude, 2026-11-07 (R2)
clear; close all; OUT = fileparts(mfilename('fullpath'));
load(fullfile(OUT, 'plot_distance_r9style_range.mat'), 'Zc', 'lat');
R1G = fullfile(fileparts(fileparts(OUT)), 'Extended_data_figures', 'FigS7_grid_cell_scatter');
mk = load(fullfile(R1G, 'topo_biome_mask_nemo.mat'), 'valid_lonlat'); mask = logical(mk.valid_lonlat)';
dn = load('I:\Prj\WOA23\nitracline\zno3_woa23_3muM_nofix.mat'); zno3 = double(dn.zno3(:, [181:360, 1:180], :));
assert(isequal(double(dn.lat(:)), lat))
NY = 20; M = nan(180, 360, 12, NY);
for y = 1:NY, d = load(['I:\Prj\MLD_1M_NEMO\one_d03\mld\mld_' num2str(2002 + y) '.mat'], 'mldG1deg'); M(:,:,:,y) = d.mldG1deg; end
mmN = [3 4 5; 6 7 8; 9 10 11; 12 1 2]; mmS = [9 10 11; 12 1 2; 3 4 5; 6 7 8]; isN = lat >= 0;
TY = nan(180, 3, 4);                                 % zonal mean: avg distance closer years, farther years, climatology check
for s = 1:4
    Ms = nan(180, 360, NY); zs = nan(180, 360);
    Ms(isN,:,:) = squeeze(mean(M(isN,:,mmN(s,:),:), 3)); Ms(~isN,:,:) = squeeze(mean(M(~isN,:,mmS(s,:),:), 3));
    zs(isN,:) = mean(zno3(isN,:,mmN(s,:)), 3, 'omitnan'); zs(~isN,:) = mean(zno3(~isN,:,mmS(s,:)), 3, 'omitnan');
    D = reshape(zs, [], 1) - reshape(Ms, [], NY);
    v = all(isfinite(D), 2) & reshape(mask, [], 1) & abs(repmat(lat, 360, 1)) <= 40;   % same cells as Zc
    cl = mean(D, 2); lo = D; lo(D >= cl) = NaN; hi = D; hi(D <= cl) = NaN;
    Q = nan(size(D,1), 3); Q(v,:) = [mean(lo(v,:), 2, 'omitnan'), mean(hi(v,:), 2, 'omitnan'), cl(v)];
    for j = 1:3, TY(:, j, s) = mean(reshape(Q(:,j), 180, 360), 2, 'omitnan'); end
end
assert(max(abs(TY(:,3,:) - Zc(:,1,:)), [], 'all') < 1e-6)   % climatology identical to 10a
seasName = {'Spring','Summer','Autumn','Winter'};
season_colors = [0.784 0.498 0.000; 0.000 0.620 0.451; 0.30 0.30 0.30; 0.000 0.447 0.698];   % as R9
lat_ticks = -40:10:40; lat_tick_labels = arrayfun(@lat_label, lat_ticks, 'UniformOutput', false);
w = abs(lat) <= 40; L = lat(w);
fig = figure('Units', 'centimeters', 'Position', [2 2 14 16], 'Color', 'w');
tl = tiledlayout(fig, 1, 1, 'Padding', 'compact'); ax = nexttile(tl); hold(ax, 'on');
for s = 1:4                                   % light bands first so lines stay on top
    c = season_colors(s,:); q = Zc(w,:,s); ty = TY(w,:,s); ok = all(isfinite(q), 2) & all(isfinite(ty), 2); q = q(ok,:); ty = ty(ok,:); Ly = L(ok);
    fill(ax, [ty(:,1); flipud(ty(:,2))], [Ly; flipud(Ly)], c, 'FaceAlpha', 0.22, 'EdgeColor', 'none');
    k = abs(Ly) > 10; fprintf('%s 10-40 N+S medians: closer years %.1f m below clim, farther years %.1f m above clim, |trend change| %.1f m\n', seasName{s}, ...
        median(q(k,1) - ty(k,1)), median(ty(k,2) - q(k,1)), median(abs(q(k,3) - q(k,2))));
end
hs = gobjects(1,4);
for s = 1:4
    c = season_colors(s,:); q = Zc(w,:,s); ok = all(isfinite(q), 2); q = q(ok,:); Ly = L(ok);
    fill(ax, [q(:,2); flipud(q(:,3))], [Ly; flipud(Ly)], c, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
    plot(ax, q(:,1), Ly, '-', 'Color', c * 0.7, 'LineWidth', 1.2); hs(s) = plot(ax, nan, nan, '-', 'Color', c, 'LineWidth', 3);
end
ylim(ax, [-40 40]); yticks(ax, lat_ticks); yticklabels(ax, lat_tick_labels); xlim(ax, [0 200]);
xlabel(ax, 'Nitracline - MLD (m)'); ylabel(ax, 'Latitude'); box(ax, 'on'); ax.FontSize = 12; ax.Toolbar.Visible = 'off';
g = [0.45 0.45 0.45];
l1 = plot(ax, nan, nan, '-', 'Color', g * 0.7, 'LineWidth', 1.2);
l2 = fill(ax, nan(3,1), nan(3,1), g, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
l3 = fill(ax, nan(3,1), nan(3,1), g, 'FaceAlpha', 0.22, 'EdgeColor', 'none');
lg = legend(ax, [hs l1 l2 l3], [seasName, {'Climatology (2003 to 2022)', 'Change due to MLD trend (2003 to 2022)', 'Change due to MLD interannual variability'}], ...
    'Box', 'off', 'FontSize', 10, 'Location', 'east');
exportgraphics(fig, fullfile(OUT, 'fig_distance_r9style_onepanel.png'), 'Resolution', 300);

function s = lat_label(v)
    if v < 0, s = sprintf('%d%cS', -v, char(176)); elseif v > 0, s = sprintf('%d%cN', v, char(176)); else, s = sprintf('0%c', char(176)); end
end
