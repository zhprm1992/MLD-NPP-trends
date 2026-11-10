% Used for: SI Fig. S5, called with 'sig005'
% Original: Manuscript\MLD\manuscript\nature\NCC\R1\code for revision\mld_npp_relationship_biome\plot_biome_sensitivity_iqr_ensemble.m (data paths updated to I:\Prj; computation unchanged)
% plot_biome_sensitivity_iqr_ensemble.m
%
% Ensemble-NPP counterpart to plot_biome_sensitivity_iqr.m. Identical
% design (median bar + IQR error bar + jittered per-pair scatter, biomes
% ordered north-to-south by true mean latitude), reading
% senMLD_biome_ensemble_<suffix>.mat instead of the CbPM file.
%
% Haipeng Zhao
% 2026-09-12

function plot_biome_sensitivity_iqr_ensemble(suffix)
OUT_DIR = fileparts(mfilename('fullpath'));
d = load(fullfile(OUT_DIR, sprintf('senMLD_biome_ensemble_%s.mat', suffix)));
senMLD = d.senMLD; senPt = d.senPt; biomes = d.biomes; season_names = d.season_names;
sig_level = d.sig_level;

NBIOME = numel(biomes);
sen_all = senMLD(:,:,1);
sen_all(isnan(sen_all)) = 0;

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
    lo = flip(sen - p25_all(ord, s));
    hi = flip(p75_all(ord, s) - sen);
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

xlabel(t, sprintf('NPP-MLD sensitivity, ensemble mean NPP (median [IQR] of significant grid-cell pairs, p < %.2f)', sig_level), 'FontWeight','bold')

out_file = fullfile(OUT_DIR, sprintf('fig_biome_sensitivity_ensemble_%s_iqr.png', suffix));
exportgraphics(fig, out_file, 'Resolution', 300);
fprintf('Saved: %s\n', out_file);
end
