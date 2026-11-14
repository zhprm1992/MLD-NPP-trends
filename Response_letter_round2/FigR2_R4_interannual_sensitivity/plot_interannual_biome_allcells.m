% Used for: round-2 letter Fig. R4 (top: Extended Data Fig. 9a samples, bottom: all grid cells)
% Original: Manuscript\MLD\manuscript\nature\NCC\R2\code for revision\interannual_biome\plot_interannual_biome_allcells.m (data paths updated to I:\Prj; computation unchanged)
% plot_interannual_biome_allcells.m (2026-11-10), internal figure for the meeting.
% Interannual GLMM biome slopes (biome as group, detrended percent anomalies), style of Figure R4 bottom row.
%   top    : Figure R4 samples (grid cell x month, both trends p < 0.05), interannual_biome_pct.mat
%   bottom : all biome grid cells and months, no significance filter, interannual_biome_allcells.mat
% Bar = biome slope, whiskers = approx 95% CI, asterisk = CI excludes 0 (no singular biome fits in Figure R4,
% bottom row not checked). No per-sample circles. READ ONLY on inputs. Haipeng Zhao / Claude, 2026-11-10.
clear; clc; OUT = fileparts(mfilename('fullpath'));
F = load(fullfile(OUT, 'interannual_biome_pct.mat'), 'Tab'); F = F.Tab;
A = load(fullfile(OUT, 'interannual_biome_allcells.mat'), 'A'); A = A.A;
seasName = {'Spring','Summer','Autumn','Winter'};
models = {'CbPM','VGPM','CAFE'}; c = [0 114 178; 0 158 115; 213 94 0] / 255;
ordB = {'NA SPSS','NP SPSS','NA STSS','NP STSS','NA STPS','NP STPS','EP EQU','A EQU','WP EQU','IND STPS','SA STPS','SP STPS','SO STSS','SO SPSS'};
NB = numel(ordB); XL = [-2.6 1];
figure('Units', 'centimeters', 'Position', [2 2 32 36], 'Color', 'w');
t = tiledlayout(2, 4); t.Padding = 'compact'; t.TileSpacing = 'compact';
for row = 1:2
    if row == 1, D = F; else, D = A; end
    for s = 1:4
        nexttile; sen = nan(NB, 3); lo = sen; hi = sen;
        for k = 1:3
            for u = 1:NB
                q = D(strcmp(D.npp, models{k}) & strcmp(D.season, seasName{s}) & strcmp(D.biome, ordB{u}), :);
                if isempty(q) || q.n == 0, continue; end
                sen(u,k) = q.slope; lo(u,k) = q.ci_lo; hi(u,k) = q.ci_hi;
            end
        end
        has = isfinite(sen); h = barh(1:NB, flip(sen, 1)); drawnow
        for k = 1:3, h(k).FaceColor = c(k,:); h(k).EdgeColor = c(k,:); h(k).LineWidth = 0.6; h(k).FaceAlpha = 0.85; h(k).BarWidth = 0.9; end
        hold on
        for k = 1:3
            ypos = flip(h(k).XData + h(k).XOffset); idx = has(:,k);
            errorbar(sen(idx,k), ypos(idx)', sen(idx,k) - lo(idx,k), hi(idx,k) - sen(idx,k), 'horizontal', 'LineStyle', 'none', 'Color', 'k', 'LineWidth', 1, 'CapSize', 3);
            addStars(sen(:,k), lo(:,k), hi(:,k), ypos, XL);
        end
        if s == 1
            if row == 2, lg = legend(h, models, 'Orientation', 'horizontal'); lg.Layout.Tile = 'south'; end
            set(gca, 'YTick', 1:NB, 'YTickLabel', flip(ordB))
        else, set(gca, 'YTick', 1:NB, 'YTickLabel', ''); end
        title(seasName{s})
        text(0, 1.005, ['(' char('a' + 4 * (row - 1) + s - 1) ')'], 'Units', 'normalized', 'VerticalAlignment', 'bottom', 'FontWeight', 'bold', 'FontSize', 11)
        xlim(XL); ylim([0.4 NB + 0.6]); xticks(-2.5:0.5:1); box on; set(gca, 'FontSize', 10, 'LineWidth', 0.8)
        xline(0, 'k-', 'LineWidth', 0.6, 'HandleVisibility', 'off')
    end
end
exportgraphics(gcf, fullfile(OUT, 'fig_interannual_biome_allcells_vs_R4.png'), 'Resolution', 200)
function addStars(s, lo, hi, y, xl)
    d = 0.04 * diff(xl);
    for i = 1:numel(s)
        if ~isfinite(s(i)) || ~(lo(i) > 0 || hi(i) < 0), continue; end
        if s(i) > 0, x = hi(i) + d; else, x = lo(i) - d; end
        x = min(max(x, xl(1) + d), xl(2) - d);
        text(x, y(i), '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
    end
end
