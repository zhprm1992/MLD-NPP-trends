% Used for: round-2 letter Fig. R2 (top: Fig. 4 samples, bottom: all latitude rows)
% Original: Manuscript\MLD\manuscript\nature\NCC\R2\code for revision\interannual_fig4\plot_glmm_allrows_vs_fig4.m (data paths updated to I:\Prj; computation unchanged)
% plot_glmm_allrows_vs_fig4.m (2026-11-10), internal figure for the meeting (Nicolas's comment on representativeness).
% Interannual GLMM band slopes (band as group, detrended percent anomalies), style of Figure R2 bottom row.
%   top    : Figure R2 samples (Fig. 4 samples, paired significant trends p < 0.1), glmm_band_group_pct.mat (main)
%   bottom : all latitude rows and months 60S-60N, no significance filter, glmm_band_group_allrows.mat
% Bar = band slope, whiskers = approx 95% CI, asterisk = CI excludes 0. No per-sample circles.
% Top row: no asterisk for singular fits, as in Figure R2. Bottom row: singular fits not checked.
% READ ONLY on inputs. Haipeng Zhao / Claude, 2026-11-10.

clear; clc;
OUT = fileparts(mfilename('fullpath'));
G = load(fullfile(OUT, 'glmm_band_group_pct.mat'), 'Tab'); G = G.Tab(strcmp(G.Tab.model, 'main'), :);
A = load(fullfile(OUT, 'glmm_band_group_allrows.mat'), 'A'); A = A.A;
SG = load(fullfile(OUT, 'glmm_singular_all_figures.mat'), 'S'); SG = SG.S(strcmp(SG.S.figure, 'fig01/03'), :);
seasName = {'Spring','Summer','Autumn','Winter'};
bandName = {'40-60N','20-40N','0-20N','0-20S','20-40S','40-60S'};
models = {'CbPM','VGPM','CAFE'};
c = [0 114 178; 0 158 115; 213 94 0] / 255;
XL = [-2 2];

figure('Units', 'centimeters', 'Position', [2 2 30 24], 'Color', 'w');
t = tiledlayout(2, 4); t.Padding = 'compact'; t.TileSpacing = 'compact';
for row = 1:2
    if row == 1, D = G; else, D = A; end
    for s = 1:4
        nexttile
        sen = nan(6, 3); lo = sen; hi = sen;
        for k = 1:3
            for u = 1:6
                q = D(strcmp(D.npp, models{k}) & strcmp(D.season, seasName{s}) & strcmp(D.band, bandName{u}), :);
                if isempty(q) || q.n_samples == 0, continue; end
                sen(u,k) = q.slope; lo(u,k) = q.ci_lo; hi(u,k) = q.ci_hi;
            end
        end
        has = isfinite(sen);
        h = barh(1:6, flip(sen, 1));
        drawnow
        for k = 1:3
            h(k).FaceColor = c(k,:); h(k).EdgeColor = c(k,:); h(k).LineWidth = 0.8; h(k).FaceAlpha = 0.85; h(k).BarWidth = 0.9;
        end
        hold on
        for k = 1:3
            ypos = flip(h(k).XData + h(k).XOffset);
            idx = has(:,k);
            errorbar(sen(idx,k), ypos(idx)', sen(idx,k) - lo(idx,k), hi(idx,k) - sen(idx,k), 'horizontal', ...
                'LineStyle', 'none', 'Color', 'k', 'LineWidth', 1, 'CapSize', 3);
            singular = row == 1 && any(SG.singular & strcmp(SG.npp, models{k}) & strcmp(SG.season, seasName{s}));
            if ~singular, addStars(sen(:,k), lo(:,k), hi(:,k), ypos, XL); end
        end
        if s == 1
            if row == 2, lg = legend(h, models, 'Orientation', 'horizontal'); lg.Layout.Tile = 'south'; end
            set(gca, 'YTickLabel', {'40 - 60S', '20 - 40S', '0 - 20S', '0 - 20N', '20 - 40N', '40 - 60N'})
        else
            set(gca, 'YTickLabel', '')
        end
        title(seasName{s})
        text(0, 1.005, ['(' char('a' + 4 * (row - 1) + s - 1) ')'], 'Units', 'normalized', ...
            'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontWeight', 'bold', 'FontSize', 11)
        xlim(XL); xticks(-2:0.5:2); box on; set(gca, 'FontSize', 10, 'LineWidth', 0.8)
        xline(0, 'k-', 'LineWidth', 0.6, 'HandleVisibility', 'off')
    end
end
exportgraphics(gcf, fullfile(OUT, 'fig_glmm_allrows_vs_fig4.png'), 'Resolution', 200)

function addStars(s, lo, hi, y, xl)
    d = 0.04 * diff(xl);
    for i = 1:numel(s)
        if ~isfinite(s(i)) || ~(lo(i) > 0 || hi(i) < 0), continue; end
        if s(i) > 0, x = hi(i) + d; else, x = lo(i) - d; end
        x = min(max(x, xl(1) + d), xl(2) - d);
        text(x, y(i), '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
    end
end
