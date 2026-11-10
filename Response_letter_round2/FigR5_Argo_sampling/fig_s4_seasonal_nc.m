% Used for: round-2 letter Fig. R5 (Argo netCDF, as in the manuscript)
% Original: Manuscript\MLD\manuscript\nature\NCC\R2\code for revision\argo_sampling\fig_s4_seasonal_nc.m (data paths updated to I:\Prj; computation unchanged)
% fig_s4_seasonal_nc.m  Figure R5 with the netCDF Argo used by the manuscript (2026-11-09); otherwise identical to fig_s4_seasonal.m.
% Figure R5 (R2 comment 3): seasonal version of manuscript Fig. S4 (Code\fs4_lat_mld_sub_fully.m),
% seasonal means as in Fig. 3. Columns: MLD-Argo, MLD-NEMO subsampled with the Argo gap mask, MLD-NEMO fully sampled.
% Statistics as in Fig. 5 (Code\LAT_vs_MLD_mkTrends_sub_full_annual_v2.m). Inputs READ ONLY. HZ / Claude, 2026-11-09.
clear; clc;
OUT = fileparts(mfilename('fullpath'));
TB = 'D:\Users\zhp19\OneDrive - Indiana State University\Research\MATLAB toolbox';
addpath(fullfile(TB, 'ktaub'), fullfile(TB, 'Colormap', 'crameri'), ...
    fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'helper_functions'));
diary(fullfile(OUT, 'fig_s4_seasonal_nc_log.txt')); diary on; fprintf('Start %s\n', datestr(now));
Y = 2003:2022; NY = 20;
A = nan(180,360,12,NY,'single'); M = A;
for y = 1:NY
    fd = dir(['I:\Prj\MLD_ARGO\one\0.03sigmaTheta_ml_' num2str(Y(y)) '.nc\*.nc']); assert(numel(fd) == 12);
    for j = 1:12, A(:,:,j,y) = ncread(fullfile(fd(j).folder, fd(j).name), 'MLD')'; end   % netCDF (lon,lat) -> (lat,lon), as the manuscript scripts
    m = load(['I:\Prj\MLD_1M_NEMO\one_d03\mld\mld_' num2str(Y(y)) '.mat'], 'mldG1deg'); M(:,:,:,y) = m.mldG1deg;
end
d = load('I:\Prj\npp_cbpm_1M_1dg_zno3mMU\grid_npp_nemo\npp_2003.mat', 'lat'); lat = double(d.lat(:));
idx = lat <= 60 & lat >= -60; latP = lat(idx);
Ms = M; Ms(isnan(A)) = NaN;
mmN = [3 4 5; 6 7 8; 9 10 11; 12 1 2]; mmS = [9 10 11; 12 1 2; 3 4 5; 6 7 8]; sN = {'Spring','Summer','Autumn','Winter'};
fig = figure('Units', 'centimeters', 'Position', [2 2 18 23], 'Color', 'w');
tl = tiledlayout(fig, 5, 3, 'TileSpacing', 'tight', 'Padding', 'compact');
cmp = crameri('lajolla', 10); colNames = {'MLD-Argo', 'MLD-NEMO subsampled', 'MLD-NEMO fully sampled'};
for s = 1:4
    Z = cell(1,3);
    for v = 1:3
        X = A; if v == 2, X = Ms; elseif v == 3, X = M; end
        z = nan(180, NY);
        for r = 1:180
            mm = mmS(s,:); if lat(r) > 0, mm = mmN(s,:); end
            z(r,:) = squeeze(mean(mean(X(r,:,mm,:), 3, 'omitnan'), 2, 'omitnan'));
        end
        Z{v} = z(idx,:);
    end
    T = cell(1,3); P = T; R = T;
    for v = 1:3, [T{v}, P{v}, R{v}] = trendRows(Z{v}); end
    for v = 1:3
        ax = nexttile(tl); axtoolbar(ax, 'Visible', 'off'); hold(ax, 'on');
        h = barh(ax, latP, T{v}, 'FaceColor', 'flat', 'EdgeColor', 'flat');
        h.CData = vals2colormap_v2(R{v}, cmp, 10, [0 1]);
        for k = find(P{v} < 0.05)', text(ax, 1.2*sign(T{v}(k)), latP(k), '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 7); end
        xlim(ax, [-1.5 1.5]); ylim(ax, [-60 60]); yticks(ax, -60:30:60); xticks(ax, -1:1:1); ax.FontSize = 7; box(ax, 'on');
        yticklabels(ax, {'60°S','30°S','0°','30°N','60°N'}); if v > 1, yticklabels(ax, []); end
        if s < 4, xticklabels(ax, []); end
        if s == 1, title(ax, colNames{v}, 'FontSize', 8, 'FontWeight', 'normal'); end
        if v == 1, text(ax, -1.4, 52, sN{s}, 'FontSize', 8); end
        if v > 1
            x = T{1}; yv = T{v}; g = isfinite(x + yv); mdl = fitlm(x(g), yv(g));
            st = [median(abs(yv(g) - x(g))), median(yv(g) - x(g)), mdl.Rsquared.Ordinary, mdl.Coefficients.Estimate(2)];
            text(ax, -1.4, -36, sprintf('vs Argo\nR^2 = %.2f\nSlope = %.2f', st(3), st(4)), 'FontSize', 7);
            fprintf('%-7s %-24s vs Argo: MAE %.2f bias %.2f R2 %.2f slope %.2f (n=%d)\n', sN{s}, colNames{v}, st, nnz(g));
        end
    end
    fprintf('%-7s significant rows (p<0.05): Argo %d, NEMO sub %d, NEMO full %d of %d\n', sN{s}, sum(P{1}<0.05), sum(P{2}<0.05), sum(P{3}<0.05), numel(latP));
end
% coverage of the Argo MLD product (same cell set as reproduce_fig5.m)
cv = squeeze(sum(sum(isfinite(A(idx,:,:,:)) & isfinite(M(idx,:,:,:)),1),2)) ./ squeeze(sum(sum(isfinite(M(idx,:,:,:)),1),2)) * 100;
axC = nexttile(tl, 13, [1 2]); axtoolbar(axC, 'Visible', 'off');
bar(axC, Y, mean(cv, 1), 0.7, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none');
ylim(axC, [0 40]); xlim(axC, [2002.4 2022.6]); box(axC, 'off'); axC.FontSize = 7;
title(axC, 'Grid cells with Argo MLD (%)', 'FontSize', 8, 'FontWeight', 'normal');
fprintf('Coverage by year (%%): %s\n', mat2str(round(mean(cv,1),1)));
axB = nexttile(tl, 15); axis(axB, 'off'); colormap(axB, cmp); caxis(axB, [0 1]);
cb = colorbar(axB, 'Location', 'southoutside'); cb.Label.String = 'R^2'; cb.FontSize = 7; cb.Ticks = 0:0.2:1;
drawnow; pB = axB.Position; cb.Position = [pB(1) + 0.1*pB(3), pB(2) + 0.45*pB(4), 0.8*pB(3), 0.08*pB(4)];
exportgraphics(fig, fullfile(OUT, 'fig_s4_seasonal_nc.png'), 'Resolution', 300);
for s = 1:4
    yrs = zeros(180,360); oc = false(180,360);
    for r = find(idx)'
        mm = mmS(s,:); if lat(r) > 0, mm = mmN(s,:); end
        yrs(r,:) = sum(squeeze(any(isfinite(A(r,:,mm,:)) & isfinite(M(r,:,mm,:)), 3)), 2)';
        oc(r,:) = all(squeeze(any(isfinite(M(r,:,mm,:)), 3)), 2)';
    end
    v = yrs(oc); fprintf('%-7s record length per cell: 20 yr %.1f%%, >=15 %.1f%%, >=10 %.1f%%, median %g yr\n', sN{s}, 100*mean(v==20), 100*mean(v>=15), 100*mean(v>=10), median(v));
end
fprintf('End %s\n', datestr(now)); diary off
function [pct, p, r2] = trendRows(Yr)
    n = size(Yr,1); pct = nan(n,1); p = nan(n,1); r2 = nan(n,1); t = (1:size(Yr,2))';
    for i = 1:n
        y = double(Yr(i,:))'; if nnz(isfinite(y)) < 10, continue; end
        [~,~,~,sig,~,~,~,sen] = ktaub([t y], 0.9999, 0);
        line = median(y, 'omitnan') + sen*(t - t(round(numel(t)/2)));
        pct(i) = sen / line(1) * 100; p(i) = sig;
        r2(i) = 1 - sum((y - line).^2, 'omitnan') / sum((y - mean(y, 'omitnan')).^2, 'omitnan');
    end
end
