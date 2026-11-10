% Used for: SI Fig. S3 (= round-1 letter Fig. R5); reads senMLD_NEMO from senMLD_ensemble.mat
% Original: Manuscript\MLD\manuscript\nature\NCC\R1\code for revision\npp_ensemble\fs3_plot_sensitivity_NPP_MLD_ensemble.m (data paths updated to I:\Prj; computation unchanged)
% Plot NPP (ensemble: VGPM+CAFE+CbPM) sensitivity to MLD (NEMO) vs Lat (median version)
% Adapted from code for revision/npp_cafe/fs3_plot_sensitivity_NPP_MLD.m
% Haipeng Zhao
% 7/12/2026
%
% SUPERSEDED 2026-07-15 by fs3_plot_sensitivity_NPP_MLD_ensemble_allModels.m
% (all 3 MLD products, one bar group per lat-band like npp_cafe's
% ..._allModels.m; this script was NEMO-only and reads a bare `senMLD`
% variable that no longer exists in senMLD_ensemble.mat). Kept for
% history; do not re-run as-is.

clear
clc

load I:\Prj\npp_ensemble_1M_1dg\senMLD_ensemble.mat
senMLD = senMLD_NEMO; senPt = senPt_nemo;  % the file now stores the NEMO results under these names (see header)
figure('Units','centimeters','Position',[2 2 30 15]);
t = tiledlayout(1, 4);
t.Padding = 'compact';

c = [120, 120, 120] / 255;

senMLD(isnan( senMLD )) = 0;

seasonNames = {'Spring', 'Summer', 'Autumn', 'Winter'};
ytl = {'40 - 60S', '20 - 40S', '0 - 20S', '0 - 20N', '20 - 40N', '40 - 60N'};

for s = 1 : 4
    nexttile
    sen = senMLD(:, s, 1);
    err = senMLD(:, s, 2);

    y = 1 : 6;
    x = flip( sen, 1 );
    h = barh( y, x );
    h.BarWidth = 0.5;
    drawnow

    h.FaceColor = c;
    h.EdgeColor = c;
    h.LineWidth = 0.8;
    h.FaceAlpha = 0.9;

    hold on;
    ypos = flip( h.XData + h.XOffset );
    for tt = 1 : 6
        pt = senPt{tt, s};
        nppPt = pt(:, 1) ./ pt(:, 2);
        nppPt = nppPt( ~isinf( nppPt ) );
        if length( nppPt ) == 1
            yposPt = ypos(tt);
        elseif isempty( nppPt )
            continue
        else
            yposPt = linspace( ypos(tt) - h.BarWidth / 18, ypos(tt) + h.BarWidth / 18, length( nppPt ) );
        end
        scatter( nppPt, yposPt, 20, 'MarkerFaceColor', c, 'MarkerFaceAlpha', 0.1, 'MarkerEdgeColor', c, 'LineWidth', 1 )
    end
    % draw errorbar
    idxErr = err ~= 0;
    errorbar( sen( idxErr ), ypos(idxErr)', err(idxErr), 'horizontal', 'LineStyle', 'none', ...
        'Color', 'k', 'LineWidth', 1 );

    title( seasonNames{s} )
    if s == 1
        set( gca, 'YTickLabel', ytl )
    else
        set( gca, 'YTickLabel', '' )
    end
    xlim( [-5 5] )
    set( gca, 'XTick', -5:2.5:5 )
end
