% Used for: Supplementary Fig. 2 (= round-1 letter Fig. R4; legend moved by hand)
% Original: Manuscript\MLD\manuscript\nature\NCC\R1\code for revision\mld_d125\f4_plot_sensitivity_NPP_MLD_d125_vs_003.m (data paths updated to I:\Prj; computation unchanged)
% f4_plot_sensitivity_NPP_MLD_d125_vs_003.m
%
% Reviewer-1 MLD-threshold sensitivity check for the manuscript's actual
% Figure 4 (NPP/MLD sensitivity vs latitude): same barh + per-point
% scatter + errorbar style as Code\f4_plot_sensitivity_NPP_MLD.m, but 2
% groups (0.03 kg/m^3, established; 0.125 kg/m^3, Reviewer-1 threshold)
% instead of 3 products (NEMO/ECCO/HYCOM) -- everything else (6 lat bands
% x 4 season panels, median sensitivity bar with std errorbar, individual
% lat-band/month ratio points overlaid) is unchanged from the original.
%
% Data: senMLD_v2_d125_vs_003.mat (senMLD2_d125_vs_003_ot.m), built via
% the identical two-stage MLD pipeline (biomes_1M_d125_ot.m ->
% monthlyLat_d125_ot.m -> compute_mld_trends_monthly_d125.m) as the
% manuscript's own 0.03 threshold, per user request 2026-07-26.
%
% Output: fig4_sensitivity_d125_vs_003.png (this folder).
%
% Haipeng Zhao
% 2026-07-26

clear; close all; clc;
DIR = fileparts(mfilename('fullpath'));
load(fullfile(DIR,'senMLD_v2_d125_vs_003.mat'));

figure('Units','centimeters','Position',[2 2 30 15]);
t = tiledlayout(1, 4);
t.Padding = 'compact';

cmp = slanCM( 'tab20c' );
cmp = unique( cmp, "rows", 'stable' );

c(1, :) = [120, 120, 120] / 255;   % 0.03 (established)
c(2, :) = [230, 159, 0]   / 255;   % 0.125 (Reviewer-1 threshold)

senMLD_NEMO_003(isnan( senMLD_NEMO_003 ))   = 0;
senMLD_NEMO_d125(isnan( senMLD_NEMO_d125 )) = 0;

season_names = {'Spring','Summer','Autumn','Winter'};
ytick_full = {'40 - 60S', '20 - 40S', '0 - 20S', '0 - 20N', '20 - 40N', '40 - 60N'};

for si = 1:4
    nexttile

    sen = cat( 2, senMLD_NEMO_003(:, si, 1), senMLD_NEMO_d125(:, si, 1) );
    err = cat( 2, senMLD_NEMO_003(:, si, 2), senMLD_NEMO_d125(:, si, 2) );

    y = 1 : 6;
    x = flip( sen, 1 );
    h = barh( y, x );
    drawnow
    for i = 1 : 2
        h(i).FaceColor = c(i, :);
        h(i).EdgeColor = c(i, :);
        h(i).LineWidth = 0.8;
        h(i).FaceAlpha = 0.9;
    end

    hold on;
    for k = 1 : size( x, 2 )
        ypos = flip( h(k).XData + h(k).XOffset );
        for tt = 1 : 6
            if k == 1
                pt = senPt_nemo_003{tt, si};
            else
                pt = senPt_nemo_d125{tt, si};
            end
            nppPt = pt(:, 1) ./ pt(:, 2);
            nppPt = nppPt( ~isinf( nppPt ) );
            if isempty(nppPt)
                continue
            end
            if length( nppPt ) == 1
                yposPt = ypos(tt);
            else
                yposPt = linspace( ypos(tt) - h(k).BarWidth / 18, ypos(tt) + h(k).BarWidth / 18, length( nppPt ) );
            end
            scatter( nppPt, yposPt, 20, 'MarkerFaceColor', cmp(k, :), 'MarkerFaceAlpha', 0.1, 'MarkerEdgeColor', c(k, :), 'LineWidth', 1 )
        end
        idx = err(:, k) ~= 0;
        errorbar( sen( idx, k ), ypos(idx)', err(idx, k), 'horizontal', 'LineStyle', 'none', ...
            'Color', 'k', 'LineWidth', 1 );
    end

    if si == 1
        h_leg = h;  % bar handles, same series order in every panel -- used for the shared legend below
        set( gca, 'YTickLabel', ytick_full )
    else
        set( gca, 'YTickLabel', '' )
    end
    xlim( [-3 3] )
    set( gca, 'XTick', -3:1:3, 'FontSize', 10, 'XTickLabelRotation', 0 )
end

lgd = legend( h_leg, {'0.03 kg m^{-3}', '0.125 kg m^{-3}'}, 'Orientation', 'horizontal' );
lgd.Layout.Tile = 'south';
lgd.Box = 'off';
lgd.FontSize = 10;

out_png = fullfile(DIR,'fig4_sensitivity_d125_vs_003.png');
exportgraphics(gcf, out_png, 'Resolution', 200);
fprintf('Figure saved: %s\n', out_png);
