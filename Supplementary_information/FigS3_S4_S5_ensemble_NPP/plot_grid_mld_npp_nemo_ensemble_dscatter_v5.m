% Used for: SI Fig. S4 (uses helper_functions/york_reg.m)
% Original: Manuscript\MLD\manuscript\nature\NCC\R1\code for revision\mld_npp_relationship_grid\plot_grid_mld_npp_nemo_ensemble_dscatter_v5.m (data paths updated to I:\Prj; computation unchanged)
% plot_grid_mld_npp_nemo_ensemble_dscatter_v5.m
%
% v3 of plot_grid_mld_npp_nemo_ensemble_dscatter.m. One correction only:
%
%   WEIGHTING CONSISTENCY. v1 area-weighted (cos-lat) only the York line
%   (via scaled input uncertainties); OLS (fitlm) and RMA (gmregress)
%   were called unweighted, so every 1x1-deg cell counted equally
%   regardless of true ocean area. v3 area-weights all three: fitlm gets
%   'Weights', and RMA is replaced by a weighted-moments version (local
%   function below), so the three lines differ only in their
%   errors-in-variables assumption, not also in spatial weighting.
%
%   Density coloring is the ORIGINAL per-panel dscatter.m, unchanged: for
%   reading density within a single panel (not comparing brightness
%   across panels of different n) per-panel normalization is exactly
%   right, so no shared-scale variant is used here (that was tried in
%   v2 and dropped as unnecessary for this purpose).
%
% Haipeng Zhao
% 2026-09-11

clear; close all; clc;
OUT_DIR = fileparts(mfilename('fullpath'));
load(fullfile(OUT_DIR,'grid_trends_mld_npp_ensemble_nemo.mat'));   % lat,lon,season_names,NPP,MLD

NLAT = numel(lat); NLON = numel(lon);
w_lat = repmat(cosd(lat), 1, NLON);

% ---- paper's open-ocean mask: 21 biomes + topo <= -1000 m -------------
% Same convention as Code/npp_cbpm/monthlyLat_nemo_ot.m (loop over all
% biome files, drop cells with topo > -1000 m), precomputed in
% mld_mechanism/topo_biome_mask_nemo.mat. v3 and earlier applied NO mask
% at all, so shelf/coastal cells that every other analysis in the paper
% excludes were being included here.
mk = load(fullfile(OUT_DIR,'..','..','Extended_data_figures','FigS7_grid_cell_scatter','topo_biome_mask_nemo.mat'), ...
          'valid_lonlat','lat');
ocean_mask = logical(mk.valid_lonlat)';   % (360 lon,180 lat) -> (180 lat,360 lon)
assert(isequal(size(ocean_mask),[NLAT NLON]), 'mask size mismatch');
assert(max(abs(double(mk.lat(:))-lat(:))) < 1e-6, 'mask latitude vector mismatch');
% v5: restrict to 60S-60N, the paper's focus domain. The lat-row pipeline
% (monthlyLat_nemo_ot.m) does NOT apply this cut, so note the difference if
% these panels are compared directly against Figure 4's latitude coverage.
ocean_mask = ocean_mask & (abs(lat(:)) <= 60);
fprintf('open-ocean mask (biome+topo, 60S-60N): %d grid cells retained; per-season retention vs cells-with-data below\n', sum(ocean_mask,'all'));

fig = figure('Units','centimeters','Position',[2 2 42 12],'Color','w');
t = tiledlayout(fig,1,4,'Padding','compact','TileSpacing','loose');

for s = 1:4
    ax = nexttile(t); hold(ax,'on');

    x  = double(MLD.trend_pct(:,:,s));  ux = double(MLD.trend_se(:,:,s));
    y  = double(NPP.trend_pct(:,:,s));  uy = double(NPP.trend_se(:,:,s));
    valid = isfinite(x) & isfinite(y) & isfinite(ux) & isfinite(uy) & ux>0 & uy>0 ...
        & ocean_mask;   % v4: paper's biome+topo open-ocean mask
    n = sum(valid,'all');
    n_hasdata = sum(isfinite(x) & isfinite(y) & isfinite(ux) & isfinite(uy) & ux>0 & uy>0, 'all');

    w = w_lat; w(~valid) = NaN;
    w = w .* sum(valid,'all') ./ sum(w,'all','omitnan');

    Xv = x(valid); Yv = y(valid); uXv = ux(valid); uYv = uy(valid); Wv = w(valid);

    % --- York (unchanged: already area-weighted via scaled uncertainties)
    [york_a, york_b, ~, york_sb] = york_reg(Xv', Yv', uXv'./sqrt(Wv'), uYv'./sqrt(Wv'));

    % --- OLS, now area-weighted (was unweighted in v1)
    mdl = fitlm(Xv, Yv, 'Weights', Wv);
    ols_a = mdl.Coefficients{1,1}; ols_b = mdl.Coefficients{2,1}; ols_sb = mdl.Coefficients{2,2};

    % --- RMA/SMA, now area-weighted (was unweighted gmregress in v1)
    [sma_a, sma_b] = local_weighted_rma(Xv, Yv, Wv);

    [rr,pp] = corr(Xv,Yv);
    fprintf('%-6s ENSEMBLE-NEMO  n=%6d  York b=%+6.3f+-%.3f  OLS(w) b=%+6.3f+-%.3f  RMA(w) b=%+6.3f  r=%+.3f (p=%.1e)\n', ...
        season_names{s}, n, york_b, york_sb, ols_b, ols_sb, sma_b, rr, pp);
    fprintf('        mask retains %d of %d cells with data (%.1f%%)\n', ...
        n, n_hasdata, 100*n/n_hasdata);

    % --- density scatter: original dscatter.m, per-panel normalized ---
    axes(ax); %#ok<LAXES>
    dscatter(Xv, Yv);

    xfit = linspace(-6, 6, 200);
    colOLS = [0.80 0.10 0.10]; colSMA = [0.45 0.20 0.70];
    hOLS = plot(ax, xfit, ols_a  + ols_b*xfit,  '-', 'Color', colOLS, 'LineWidth', 1.6);
    hSMA = plot(ax, xfit, sma_a  + sma_b*xfit,  '-', 'Color', colSMA, 'LineWidth', 1.6);
    hYork = plot(ax, xfit, york_a + york_b*xfit, 'k-', 'LineWidth', 2.2);

    hPoint = plot(ax, nan, nan, 'o', 'MarkerFaceColor',[0.5 0.5 0.5], ...
        'MarkerEdgeColor','none', 'MarkerSize', 5);

    yline(ax,0,'Color',[.6 .6 .6],'LineStyle',':'); xline(ax,0,'Color',[.6 .6 .6],'LineStyle',':');
    xlim(ax,[-6 6]); ylim(ax,[-12 12]);
    xlabel(ax,'MLD trend (% yr^{-1})','FontSize',15); if s==1, ylabel(ax,'NPP trend (% yr^{-1})','FontSize',15); end
    title(ax, sprintf('%s (n=%d)', season_names{s}, n), 'FontWeight','normal','FontSize',16);
    set(ax,'FontSize',14,'Box','on');
    if s == 4
        cb = colorbar(ax); cb.Label.String = 'Density'; cb.Label.FontSize = 14; cb.FontSize = 13;
    end

    legend(ax, [hPoint hYork hOLS hSMA], ...
        {'1\times1^{\circ} grid cell', ...
         sprintf('York: %+.2f', york_b), ...
         sprintf('OLS: %+.2f', ols_b), ...
         sprintf('RMA: %+.2f', sma_b)}, ...
        'Location','northeast', 'FontSize',13, 'Box','on');
end

out_png = fullfile(OUT_DIR,'fig_grid_mld_npp_trend_regression_ENSEMBLE_NEMO_v5.png');
exportgraphics(fig, out_png, 'Resolution',200);
fprintf('\nFigure saved: %s\n', out_png);


%% ---------- local functions ----------
function [a,b] = local_weighted_rma(X,Y,W)
% Weighted geometric-mean / reduced-major-axis (SMA) regression.
% Standard unweighted RMA: b = sign(r)*std(Y)/std(X). This is the
% weighted-moments generalization, using the same area weights as OLS
% and (via scaled uncertainties) York, so all three lines are on equal
% spatial footing.
W = W(:); X = X(:); Y = Y(:);
sw = sum(W);
mx = sum(W.*X)/sw; my = sum(W.*Y)/sw;
vx = sum(W.*(X-mx).^2)/sw;
vy = sum(W.*(Y-my).^2)/sw;
cxy = sum(W.*(X-mx).*(Y-my))/sw;
r = cxy / sqrt(vx*vy);
b = sign(r) * sqrt(vy/vx);
a = my - b*mx;
end
