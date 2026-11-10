% Used for: Extended Data Fig. 9a-d, called as plot_biome_sensitivity_iqr('sig005'); panel (e) biome map added by hand
% Original: Manuscript\MLD\manuscript\nature\NCC\R1\code for revision\mld_npp_relationship_biome\plot_biome_sensitivity_iqr.m (data paths updated to I:\Prj; computation unchanged)
% plot_biome_sensitivity_iqr.m
%
% Same as plot_biome_sensitivity.m, but replaces the std-based error bar
% with an IQR-based one (25th-75th percentile of the significant-pair
% sensitivity ratios), since sensitivity is a ratio (NPP_trend/MLD_trend)
% and a handful of near-zero-MLD-trend denominators can blow up the raw
% std by orders of magnitude even though the median stays sensible
% (e.g. SO SPSS Spring: median=-1.53, std=392, IQR=[-2.52,-1.15] --
% confirmed empirically 2026-07-30). IQR pairs naturally with a median
% point estimate since both are robust (percentile-based) statistics.
%
% Haipeng Zhao
% 2026-07-30

function plot_biome_sensitivity_iqr(suffix)
OUT_DIR = fileparts(mfilename('fullpath'));
d = load(fullfile(OUT_DIR, sprintf('senMLD_biome_nemo_%s.mat', suffix)));
senMLD = d.senMLD; senPt = d.senPt; biomes = d.biomes; season_names = d.season_names;
sig_level = d.sig_level;

NBIOME = numel(biomes);
sen_all = senMLD(:,:,1);
sen_all(isnan(sen_all)) = 0;

% Compute IQR (25th/75th percentile) of the significant-pair ratios
p25_all = nan(NBIOME,4); p75_all = nan(NBIOME,4);
for bi = 1:NBIOME
    for s = 1:4
        pt = senPt{bi,s};
        nppPt = pt(:,1) ./ pt(:,2);
        nppPt = nppPt(~isinf(nppPt) & ~isnan(nppPt));
        if isempty(nppPt), continue; end
        p25_all(bi,s) = prctile(nppPt, 25);
        p75_all(bi,s) = prctile(nppPt, 75);
    end
end

% Ordered strictly north-to-south by each biome's actual mean latitude
% (computed from the real idx data 2026-07-31, not just grouped by
% regime name -- SA/SP STPS are Southern Hemisphere and belong after the
% equatorial biomes, not alongside NA/NP STPS).
biome_order = {'NA SPSS','NP SPSS','NA STSS','NP STSS','NA STPS','NP STPS', ...
    'EP EQU','A EQU','WP EQU', ...
    'IND STPS','SA STPS','SP STPS','SO STSS','SO SPSS'};
[~, ord] = ismember(biome_order, biomes);
assert(all(ord > 0), 'biome name mismatch');

fig = figure('Units','centimeters','Position',[2 2 34 18]);
t = tiledlayout(1, 4);
t.Padding = 'compact';

c_bar = [204, 121, 167] / 255;
c_pt  = [120, 120, 120] / 255;

for s = 1:4
    nexttile
    sen = sen_all(ord, s);
    lo = flip(sen - p25_all(ord, s));   % lower error (median - p25)
    hi = flip(p75_all(ord, s) - sen);   % upper error (p75 - median)
    y = 1:NBIOME;
    x = flip(sen);
    h = barh(y, x, 'FaceColor', c_bar, 'EdgeColor', c_bar, 'FaceAlpha', 0.9);
    drawnow
    hold on

    for t2 = 1:NBIOME
        bidx = ord(NBIOME - t2 + 1);
        pt = senPt{bidx, s};
        nppPt = pt(:,1) ./ pt(:,2);
        nppPt = nppPt(~isinf(nppPt) & ~isnan(nppPt));
        if isempty(nppPt), continue; end
        ypos = t2 * ones(size(nppPt)) + linspace(-0.15, 0.15, numel(nppPt))';
        scatter(nppPt, ypos, 14, 'MarkerFaceColor', c_pt, 'MarkerFaceAlpha', 0.15, ...
            'MarkerEdgeColor', c_pt, 'LineWidth', 0.5)
    end

    idxE = ~isnan(lo) & ~isnan(hi);
    errorbar(x(idxE), y(idxE)', lo(idxE), hi(idxE), 'horizontal', 'LineStyle','none', 'Color','k', 'LineWidth', 1)

    title(season_names{s})
    if s == 1
        set(gca, 'YTickLabel', flip(biome_order), 'YTick', 1:NBIOME, 'FontSize', 9)
    else
        set(gca, 'YTickLabel', '')
    end
    ylim([0.5 NBIOME+0.5])
    xlim([-3 3])
    xline(0, 'Color', [0.5 0.5 0.5], 'LineStyle', ':')
end

xlabel(t, sprintf('NPP-MLD sensitivity (median [IQR] of significant grid-cell pairs, p < %.2f)', sig_level), 'FontWeight','bold')

out_file = fullfile(OUT_DIR, sprintf('fig_biome_sensitivity_nemo_%s_iqr.png', suffix));
exportgraphics(fig, out_file, 'Resolution', 300);
fprintf('Saved: %s\n', out_file);
end
