% Used for: Supplementary Fig. 1 (= round-1 letter Fig. R3)
% Original: Manuscript\MLD\manuscript\nature\NCC\R1\code for revision\mld_d125\plot_mldLat_trend_comparison_d125_vs_003_seasonal.m (data paths updated to I:\Prj; computation unchanged)
% plot_mldLat_trend_comparison_d125_vs_003_seasonal.m
%
% 1x4 panel (one per season) redesign of
% fig_mldLat_trend_comparison_d125_vs_003.png, styled after the
% manuscript's actual Figure 5 (Code\scatterMLD_ARGOvsMODEL_trends_sp.m
% and its su/au/wi siblings): scatter colored by latitude (jet colormap,
% caxis [-60 60]), 1:1 dashed reference line in the same maroon
% [0.6350 0.0780 0.1840], same [-1.5, 1.5] axis range convention. Each
% panel title is just the season name. Left-corner annotation gives R^2
% and slope (OLS fit of 0.125-threshold trend vs 0.03-threshold trend,
% via fitlm) plus the sign-agreement count for that season, per user
% request 2026-07-26.
%
% Data: mldLat_trend_comparison_d125_vs_003.csv (from
% compare_mldLat_trends_d125_vs_d03.m -- Sen's-slope/Mann-Kendall
% latitudinal MLD trend, 0.03 vs 0.125 kg/m^3 threshold, same ktaub
% methodology as the manuscript's own LAT_vs_MLDmkTrends_season.m).
%
% Output: fig_mldLat_trend_comparison_d125_vs_003_seasonal.png (this folder).
%
% Haipeng Zhao
% 2026-07-26

clear; close all; clc;
DIR = fileparts(mfilename('fullpath'));
T = readtable(fullfile(DIR,'mldLat_trend_comparison_d125_vs_003.csv'));

seasons = {'Spring','Summer','Autumn','Winter'};
lims = [-1.5 1.5];

fig = figure('Units','centimeters','Position',[2 2 34 12],'Color','w');
t = tiledlayout(fig, 1, 4, 'Padding','loose', 'TileSpacing','loose');

for si = 1:4
    ax = nexttile(t);
    m = strcmp(T.season, seasons{si});
    x = T.trend003(m); y = T.trend125(m); latv = T.lat(m);

    scatter(ax, x, y, 28, latv, 'filled');
    cmp = jet(128);
    colormap(ax, cmp);
    caxis(ax, [-60 60]);
    if si == 4
        cb = colorbar(ax);
        cb.Ticks = -60:20:60;
        cb.Label.String = 'Latitude (\circ)';
    end

    hold(ax, 'on');
    xline(ax, 0, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
    yline(ax, 0, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
    plot(ax, lims, lims, '--', 'LineWidth', 1.5, 'Color', [0.6350 0.0780 0.1840]);

    mdl = fitlm(x, y);
    intercept = mdl.Coefficients.Estimate(1);
    slope = mdl.Coefficients.Estimate(2);
    r2 = mdl.Rsquared.Ordinary;
    plot(ax, lims, intercept + slope*lims, '-', 'LineWidth', 1.5, 'Color', 'k');

    xlim(ax, lims); ylim(ax, lims);
    axis(ax, 'square');
    title(ax, seasons{si}, 'FontWeight','normal');

    n_agree = sum(T.sign_agree(m));
    n_tot = sum(m);

    txt = sprintf('R^2 = %.2f\nslope = %.2f\nn_{agree} = %d/%d', r2, slope, n_agree, n_tot);
    text(ax, lims(1)+0.05*range(lims), lims(2)-0.05*range(lims), txt, ...
        'VerticalAlignment','top', 'HorizontalAlignment','left', 'FontSize', 9);

    set(ax, 'FontSize', 10, 'Box','on');
end

xlabel(t, 'Slope of latitudinal trends of MLD, 0.03 kg m^{-3} threshold (% yr^{-1})', 'FontSize', 11);
ylabel(t, {'Slope of latitudinal trends of MLD,', '0.125 kg m^{-3} threshold (% yr^{-1})'}, 'FontSize', 11);

out_png = fullfile(DIR,'fig_mldLat_trend_comparison_d125_vs_003_seasonal.png');
exportgraphics(fig, out_png, 'Resolution', 200);
fprintf('Figure saved: %s\n', out_png);
