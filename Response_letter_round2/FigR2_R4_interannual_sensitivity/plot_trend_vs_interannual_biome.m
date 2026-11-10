% Used for: round-2 letter Fig. R4
% Original: Manuscript\MLD\manuscript\nature\NCC\R2\code for revision\interannual_biome\plot_trend_vs_interannual_biome.m (data paths updated to I:\Prj; computation unchanged)
% plot_trend_vs_interannual_biome.m  (2026-11-08, HZ: biome version of figure 03b)
% Same S8a samples as interannual_biome.m (grid cell x month, both NPP and MLD-NEMO trends p < 0.05, 14 biomes):
%   top    : trend-to-trend sensitivity per sample = NPP trend / MLD trend (% per yr), bar = median, whiskers = IQR (as S8a)
%   bottom : interannual GLMM (biome as group, detrended percent anomalies), bar = biome slope, whiskers = approx 95% CI,
%            circles = detrended OLS slope of each sample. No singular fit in the biome GLMMs (glmm_singular_all_figures).
% Both rows x axis -4 to 4 (as 03b). Check: CbPM top-row medians vs R1 senMLD_biome_nemo_sig005.mat.
% READ ONLY on inputs. Haipeng Zhao / Claude, 2026-11-08 (R2)
clear; clc; OUT = fileparts(mfilename('fullpath'));
R1B = fullfile(fileparts(fileparts(OUT)), 'Extended_data_figures', 'FigS8_biome_sensitivity');
load(fullfile(OUT, 'interannual_biome_pct.mat'), 'Tab', 'SAMP', 'biomes', 'seasName');
Tr = load(fullfile(R1B, 'grid_monthly_trends_mld_npp_nemo.mat'), 'npp_trend_pct', 'mld_trend_pct');
Tv = load(fullfile(OUT, 'grid_monthly_trends_npp_vgpm_cafe.mat'), 'vgpm_trend_pct', 'cafe_trend_pct');
S8 = load(fullfile(R1B, 'senMLD_biome_nemo_sig005.mat'), 'senMLD');
tn.CbPM = Tr.npp_trend_pct; tn.VGPM = Tv.vgpm_trend_pct; tn.CAFE = Tv.cafe_trend_pct;
models = {'CbPM','VGPM','CAFE'}; c = [0 114 178; 0 158 115; 213 94 0] / 255;
ordB = {'NA SPSS','NP SPSS','NA STSS','NP STSS','NA STPS','NP STPS','EP EQU','A EQU','WP EQU','IND STPS','SA STPS','SP STPS','SO STSS','SO SPSS'};
NB = numel(ordB); [~, bIdx] = ismember(ordB, biomes); XL = [-4 4];
% trend-to-trend sensitivity per sample (same samples as the interannual slopes)
TS = struct(); mxd = 0;
for k = 1:3
    for s = 1:4
        I = SAMP.(models{k}){s}; v = nan(size(I,1), 1);
        for i = 1:size(I,1)
            [ia, jo] = ind2sub([180 360], I(i,2));
            v(i) = double(tn.(models{k})(ia, jo, I(i,3))) / double(Tr.mld_trend_pct(ia, jo, I(i,3)));
        end
        TS.(models{k}){s} = v;
        if k == 1
            for b = 1:numel(biomes)
                w = I(:,1) == b; if ~any(w), continue; end
                d = abs(median(v(w)) - S8.senMLD(b, s, 1)); mxd = max(mxd, d);
                if d > 1e-3, fprintf('  CbPM %s %s: median %.3f vs S8a %.3f (n %d vs %d)\n', biomes{b}, seasName{s}, median(v(w)), S8.senMLD(b,s,1), sum(w), S8.senMLD(b,s,3)); end
            end
        end
    end
end
fprintf('Ground truth: max |CbPM median - S8a| = %.2g\n', mxd);
rng(1); figure('Units', 'centimeters', 'Position', [2 2 32 36], 'Color', 'w');
t = tiledlayout(2, 4); t.Padding = 'compact'; t.TileSpacing = 'compact'; nClip = [0 0];
for row = 1:2
    for s = 1:4
        nexttile; sen = nan(NB, 3); lo = sen; hi = sen; pts = cell(NB, 3);
        for k = 1:3
            I = SAMP.(models{k}){s};
            for u = 1:NB
                w = I(:,1) == bIdx(u); if ~any(w), continue; end
                if row == 1
                    v = TS.(models{k}){s}(w); sen(u,k) = median(v); lo(u,k) = prctile(v, 25); hi(u,k) = prctile(v, 75);
                else
                    q = Tab(strcmp(Tab.npp, models{k}) & strcmp(Tab.season, seasName{s}) & strcmp(Tab.biome, ordB{u}), :);
                    sen(u,k) = q.slope; lo(u,k) = q.ci_lo; hi(u,k) = q.ci_hi; v = I(w, 4);
                end
                pts{u,k} = v;
            end
        end
        has = isfinite(sen); h = barh(1:NB, flip(sen, 1)); drawnow
        for k = 1:3, h(k).FaceColor = c(k,:); h(k).EdgeColor = c(k,:); h(k).LineWidth = 0.6; h(k).FaceAlpha = 0.85; h(k).BarWidth = 0.9; end
        hold on
        for k = 1:3
            ypos = flip(h(k).XData + h(k).XOffset);
            for u = 1:NB
                p = pts{u,k}; if isempty(p), continue; end
                nClip(row) = nClip(row) + sum(abs(p) > XL(2)); p = p(abs(p) <= XL(2)); if isempty(p), continue; end
                scatter(p, ypos(u) + (rand(numel(p), 1) - 0.5) * 0.2, 6, 'MarkerFaceColor', 'w', 'MarkerFaceAlpha', 0.5, ...
                    'MarkerEdgeColor', c(k,:) * 0.55, 'MarkerEdgeAlpha', 0.5, 'LineWidth', 0.4)
            end
            idx = has(:,k);
            errorbar(sen(idx,k), ypos(idx)', sen(idx,k) - lo(idx,k), hi(idx,k) - sen(idx,k), 'horizontal', 'LineStyle', 'none', 'Color', 'k', 'LineWidth', 1, 'CapSize', 3);
            if row == 2, addStars(sen(:,k), lo(:,k), hi(:,k), ypos, XL); end
        end
        if s == 1
            if row == 1, lg = legend(h, models, 'Orientation', 'horizontal'); lg.Layout.Tile = 'south'; end
            set(gca, 'YTick', 1:NB, 'YTickLabel', flip(ordB))
        else, set(gca, 'YTick', 1:NB, 'YTickLabel', ''); end
        title(seasName{s})
        text(0, 1.005, ['(' char('a' + 4 * (row - 1) + s - 1) ')'], 'Units', 'normalized', 'VerticalAlignment', 'bottom', 'FontWeight', 'bold', 'FontSize', 11)
        xlim(XL); ylim([0.4 NB + 0.6]); xticks(-4:1:4); box on; set(gca, 'FontSize', 10, 'LineWidth', 0.8)
        xline(0, 'k-', 'LineWidth', 0.6, 'HandleVisibility', 'off')
    end
end
exportgraphics(gcf, fullfile(OUT, 'fig_trend_vs_interannual_biome.png'), 'Resolution', 300)
fprintf('Sample values beyond +-4 not drawn: top %d, bottom %d\n', nClip);
for k = 1:3
    nt = 0; nn = 0; ni = 0; nb = 0;
    for s = 1:4
        I = SAMP.(models{k}){s};
        for u = 1:NB
            w = I(:,1) == bIdx(u); if ~any(w), continue; end
            q = Tab(strcmp(Tab.npp, models{k}) & strcmp(Tab.season, seasName{s}) & strcmp(Tab.biome, ordB{u}), :);
            mt = median(TS.(models{k}){s}(w)); nb = nb + 1; nt = nt + (mt < 0); ni = ni + (q.slope < 0);
            nn = nn + (abs(q.slope) < abs(mt));
        end
    end
    fprintf('%s: biome-seasons %d | trend sens median < 0: %d | interannual GLMM < 0: %d | |interannual| < |trend|: %d\n', models{k}, nb, nt, ni, nn);
end
function addStars(s, lo, hi, y, xl)
    d = 0.04 * diff(xl);
    for i = 1:numel(s)
        if ~isfinite(s(i)) || ~(lo(i) > 0 || hi(i) < 0), continue; end
        if s(i) > 0, x = hi(i) + d; else, x = lo(i) - d; end
        x = min(max(x, xl(1) + d), xl(2) - d);
        text(x, y(i), '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
    end
end
