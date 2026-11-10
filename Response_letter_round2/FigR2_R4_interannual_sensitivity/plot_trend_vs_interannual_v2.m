% Used for: round-2 letter Fig. R2
% Original: Manuscript\MLD\manuscript\nature\NCC\R2\code for revision\interannual_fig4\plot_trend_vs_interannual_v2.m (data paths updated to I:\Prj; computation unchanged)
% plot_trend_vs_interannual_v2.m  (2026-11-08) as plot_trend_vs_interannual.m, but no asterisk for slopes from singular GLMM fits (glmm_singular_all_figures.mat).
% Two-row comparison on the same Fig. 4 samples (lat row x month, paired significant trends p < 0.1):
%   top    : trend-to-trend sensitivity (NPP trend / MLD trend, % per yr), bar = median, whiskers = IQR,
%            circles = each sample (beyond the axis not drawn, as in S8a)
%   bottom : interannual GLMM (band as group, detrended percent anomalies), bar = band slope, whiskers = approx 95% CI,
%            circles = detrended OLS slope of each sample (same as fig_glmm_band_group_pct.png)
% Inputs: trend_vs_interannual.mat, glmm_band_group_pct.mat. READ ONLY. Haipeng Zhao / Claude, 2026-11-06 (R2)

clear; clc;
OUT = fileparts(mfilename('fullpath'));
load(fullfile(OUT, 'trend_vs_interannual.mat'), 'ALL');
G = load(fullfile(OUT, 'glmm_band_group_pct.mat'), 'Tab'); G = G.Tab(strcmp(G.Tab.model, 'main'), :);
SG = load(fullfile(OUT, 'glmm_singular_all_figures.mat'), 'S'); SG = SG.S(strcmp(SG.S.figure, 'fig01/03'), :);
seasName = {'Spring','Summer','Autumn','Winter'};
bandName = {'40-60N','20-40N','0-20N','0-20S','20-40S','40-60S'};
models = {'CbPM','VGPM','CAFE'};
c = [0 114 178; 0 158 115; 213 94 0] / 255;           % NPP-model colors (not the Fig. 4 MLD colors)
XL = {[-4 4], [-4 4]};                               % HZ 2026-11-08: same x range in both rows
xt = {-4:1:4, -4:1:4};
ylab = {'Trend-to-trend sensitivity', 'Interannual sensitivity (GLMM)'};

rng(1);
figure('Units', 'centimeters', 'Position', [2 2 30 26], 'Color', 'w');
t = tiledlayout(2, 4); t.Padding = 'compact'; t.TileSpacing = 'compact';
nClip = 0;
for row = 1:2
    for s = 1:4
        nexttile
        sen = nan(6, 3); lo = sen; hi = sen; pts = cell(6, 3);
        for k = 1:3
            A = ALL.(models{k});
            for u = 1:6
                sel = A(:,1) == u & A(:,2) == s;
                if ~any(sel), continue; end
                if row == 1
                    v = A(sel, 3); sen(u,k) = median(v); lo(u,k) = prctile(v, 25); hi(u,k) = prctile(v, 75);
                else
                    q = G(strcmp(G.npp, models{k}) & strcmp(G.season, seasName{s}) & strcmp(G.band, bandName{u}), :);
                    if isempty(q) || q.n_samples == 0, continue; end
                    sen(u,k) = q.slope; lo(u,k) = q.ci_lo; hi(u,k) = q.ci_hi; v = A(sel, 4);
                end
                pts{u,k} = v;
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
            for u = 1:6
                p = pts{u,k}; if isempty(p), continue; end
                nClip = nClip + sum(abs(p) > XL{row}(2)); p = p(abs(p) <= XL{row}(2));
                if isempty(p), continue; end
                yp = ypos(u) + (rand(numel(p), 1) - 0.5) * 0.12;
                scatter(p, yp, 10, 'MarkerFaceColor', 'w', 'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', c(k,:) * 0.55, 'MarkerEdgeAlpha', 0.8, 'LineWidth', 0.5)
            end
            idx = has(:,k);
            errorbar(sen(idx,k), ypos(idx)', sen(idx,k) - lo(idx,k), hi(idx,k) - sen(idx,k), 'horizontal', ...
                'LineStyle', 'none', 'Color', 'k', 'LineWidth', 1, 'CapSize', 3);
            if row == 2 && ~any(SG.singular & strcmp(SG.npp, models{k}) & strcmp(SG.season, seasName{s})), addStars(sen(:,k), lo(:,k), hi(:,k), ypos, XL{row}); end
        end
        if s == 1
            if row == 1, lg = legend(h, models, 'Orientation', 'horizontal'); lg.Layout.Tile = 'south'; end
            set(gca, 'YTickLabel', {'40 - 60S', '20 - 40S', '0 - 20S', '0 - 20N', '20 - 40N', '40 - 60N'})
        else
            set(gca, 'YTickLabel', '')
        end
        title(seasName{s})
        text(0, 1.005, ['(' char('a' + 4 * (row - 1) + s - 1) ')'], 'Units', 'normalized', ...
            'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontWeight', 'bold', 'FontSize', 11)   % panel letters (a)-(h)
        xlim(XL{row}); xticks(xt{row}); box on; set(gca, 'FontSize', 10, 'LineWidth', 0.8)
        xline(0, 'k-', 'LineWidth', 0.6, 'HandleVisibility', 'off')
    end
end
exportgraphics(gcf, fullfile(OUT, 'fig_trend_vs_interannual_v2.png'), 'Resolution', 300)
fprintf('%d sample values beyond the axis not drawn\n', nClip);

function addStars(s, lo, hi, y, xl)
    % asterisk beyond the whisker tip (away from zero) where the approximate 95% CI excludes zero
    d = 0.04 * diff(xl);
    for i = 1:numel(s)
        if ~isfinite(s(i)) || ~(lo(i) > 0 || hi(i) < 0), continue; end
        if s(i) > 0, x = hi(i) + d; else, x = lo(i) - d; end
        x = min(max(x, xl(1) + d), xl(2) - d);
        text(x, y(i), '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
    end
end
