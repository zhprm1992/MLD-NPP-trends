% Used for: Fig. S7 (needs dscatter.m from MATLAB File Exchange)
% Original: Manuscript\MLD\manuscript\nature\NCC\R1\code for revision\mld_npp_relationship_grid\plot_grid_mld_npp_nemo_dscatter_v5.m (data paths updated to I:\Prj; computation unchanged)
% plot_grid_mld_npp_nemo_dscatter_v5.m   [MAIN TEXT FIGURE]
%
% v3 changes, relative to the original v1 (2026-07-17):
%
%   1. STALE DATA FIXED. v1's PNG (2026-07-22) predates the correction of
%      the MLD source directory from I:\Prj\MLD_1M_NEMO\one\ to
%      \one_d03\mld\ (compute_grid_seasonal_trends.m, 2026-07-25), so it
%      was drawn from superseded MLD. v3 is rerun against
%      grid_trends_mld_npp.mat as currently saved. Same staleness was
%      found and fixed in the ensemble-NEMO sibling figure.
%
%   2. ALL THREE LINES NOW AREA-WEIGHTED. v1 area-weighted only York (via
%      cos-lat-scaled input uncertainties) and left OLS/RMA unweighted.
%      That asymmetry was deliberate in v1 -- see the original header
%      below, which documents OLS/RMA as biased references illustrating
%      regression-dilution attenuation. SUPERSEDED 2026-09-11 by explicit
%      decision: OLS now takes 'Weights', and RMA is replaced by a
%      weighted-moments version (local function below), so all three
%      estimators share identical cos-lat spatial weighting and differ
%      ONLY in their errors-in-variables assumption. This also makes the
%      main-text figure consistent with the supplementary ensemble
%      figure. The attenuation point of the original design still reads
%      the same way (weighted OLS remains visibly shallower than York);
%      it is simply no longer confounded with a weighting difference.
%
%   3. Larger fonts throughout.
%
% Haipeng Zhao / 2026-09-11
%
% ---- original v1 header follows (note item 2 supersedes its OLS/RMA
%      unweighted rationale) ----
%
% Revised fig_grid_mld_npp_trend_regression.png: NEMO only (dropping
% ECCO/HYCOM from this figure for now), full point cloud (not the
% 4000-point subsample the multi-product version used) rendered as a
% density-colored scatter via dscatter.m (Eilers & Goeman 2004 smoothed-
% density method; Robert Henson's implementation, MATLAB Bioinformatics
% Toolbox / File Exchange #8430) instead of alpha-transparency, per user
% request.
%
% Three regression lines per panel, matching the comparison style already
% used in I:\Prj3\Cael_2026\model2_regression_comparison.m:
%   - York  (Model II, weighted by per-cell uncertainty + cos-lat area) --
%     PRIMARY method. Both axes have real measurement error here (each
%     cell's MLD/NPP trend has its own SE from the per-cell OLS fit), so
%     Model II is the methodologically correct choice.
%   - OLS   (Model I, unweighted) -- NOT a competing correct estimate.
%     Shown as a biased reference to demonstrate the regression-dilution
%     attenuation that results from (incorrectly) ignoring x-axis error --
%     visually justifies the York choice rather than undermining it.
%   - SMA   (Standardized Major Axis, unweighted Model II -- also called
%     Reduced Major Axis/RMA or geometric mean regression; computed via
%     gmregress.m). A genuinely different assumption from York (ignores
%     the per-cell uncertainty weighting entirely, uses only the raw
%     variance ratio) -- included as an independent Model II cross-check.
%   Deming regression considered and deliberately excluded: with a single
%   dataset-wide variance ratio it would be a cruder approximation of what
%   York already does more precisely using each cell's own uncertainty --
%   not a distinct alternative worth the extra line.
%
% Haipeng Zhao
% 2026-07-17

clear; close all; clc;
OUT_DIR = fileparts(mfilename('fullpath'));
load(fullfile(OUT_DIR,'grid_trends_mld_npp.mat'));   % lat,lon,season_names,product_names,NPP,MLD

p = find(strcmp(product_names,'nemo'));
NLAT = numel(lat); NLON = numel(lon);
w_lat = repmat(cosd(lat), 1, NLON);

% ---- paper's open-ocean mask: 21 biomes + topo <= -1000 m -------------
% Same convention as Code/npp_cbpm/monthlyLat_nemo_ot.m (loop over all
% biome files, drop cells with topo > -1000 m), precomputed in
% mld_mechanism/topo_biome_mask_nemo.mat. v3 and earlier applied NO mask
% at all, so shelf/coastal cells that every other analysis in the paper
% excludes were being included here.
mk = load(fullfile(OUT_DIR,'topo_biome_mask_nemo.mat'), ...
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

    x  = double(MLD.trend_pct(:,:,s,p));  ux = double(MLD.trend_se(:,:,s,p));
    y  = double(NPP.trend_pct(:,:,s,p));  uy = double(NPP.trend_se(:,:,s,p));

    valid = isfinite(x) & isfinite(y) & isfinite(ux) & isfinite(uy) & ux>0 & uy>0 ...
        & ocean_mask;   % v4: paper's biome+topo open-ocean mask
    n = sum(valid,'all');
    n_hasdata = sum(isfinite(x) & isfinite(y) & isfinite(ux) & isfinite(uy) & ux>0 & uy>0, 'all');

    w = w_lat; w(~valid) = NaN;
    w = w .* sum(valid,'all') ./ sum(w,'all','omitnan');

    Xv = x(valid); Yv = y(valid); uXv = ux(valid); uYv = uy(valid); Wv = w(valid);

    [york_a, york_b, ~, york_sb] = york_reg(Xv', Yv', uXv'./sqrt(Wv'), uYv'./sqrt(Wv'));

    % OLS, area-weighted (v3: was unweighted in v1/v2)
    mdl = fitlm(Xv, Yv, 'Weights', Wv);
    ols_a = mdl.Coefficients{1,1}; ols_b = mdl.Coefficients{2,1}; ols_sb = mdl.Coefficients{2,2};

    % RMA/SMA, area-weighted (v3: was unweighted gmregress in v1/v2)
    [sma_a, sma_b] = local_weighted_rma(Xv, Yv, Wv);

    [rr,pp] = corr(Xv,Yv);
    fprintf('%-6s NEMO  n=%6d  York b=%+6.3f+-%.3f  OLS(w) b=%+6.3f+-%.3f  RMA(w) b=%+6.3f  r=%+.3f (p=%.1e)\n', ...
        season_names{s}, n, york_b, york_sb, ols_b, ols_sb, sma_b, rr, pp);
    fprintf('        mask retains %d of %d cells with data (%.1f%%)\n', ...
        n, n_hasdata, 100*n/n_hasdata);

    axes(ax); %#ok<LAXES> % dscatter.m plots to gca internally (no axis-handle argument in its signature)
    dscatter(Xv, Yv);
    xfit = linspace(-6, 6, 200); % span the fixed axis width, not just this season's data range
    % All solid, differentiated by color (+ thickness for York, the primary
    % method). Tried dashed/dash-dot + color first (kept redundant line-style
    % encoding for grayscale/colorblind safety) but the user preferred the
    % cleaner all-solid look after an A/B comparison -- accessibility
    % trade-off noted and accepted deliberately, not an oversight.
    colOLS = [0.80 0.10 0.10]; colSMA = [0.45 0.20 0.70];
    hOLS = plot(ax, xfit, ols_a  + ols_b*xfit,  '-', 'Color', colOLS, 'LineWidth', 1.6);
    hSMA = plot(ax, xfit, sma_a  + sma_b*xfit,  '-', 'Color', colSMA, 'LineWidth', 1.6);
    hYork = plot(ax, xfit, york_a + york_b*xfit, 'k-', 'LineWidth', 2.2);

    % Proxy marker (not real data, NaN/NaN so nothing is actually drawn) purely
    % to give the scatter cloud its own legend entry alongside the three lines.
    hPoint = plot(ax, nan, nan, 'o', 'MarkerFaceColor',[0.5 0.5 0.5], ...
        'MarkerEdgeColor','none', 'MarkerSize', 5);

    yline(ax,0,'Color',[.6 .6 .6],'LineStyle',':'); xline(ax,0,'Color',[.6 .6 .6],'LineStyle',':');
    xlim(ax,[-6 6]); ylim(ax,[-12 12]);
    xlabel(ax,'MLD trend (% yr^{-1})','FontSize',15); if s==1, ylabel(ax,'NPP trend (% yr^{-1})','FontSize',15); end
    title(ax, sprintf('%s (n=%d)', season_names{s}, n), 'FontWeight','normal','FontSize',16);
    set(ax,'FontSize',14,'Box','on');
    if s == 4  % all 4 panels share the same [0,1] density colormap/scale -- one shared colorbar suffices
        cb = colorbar(ax); cb.Label.String = 'Density'; cb.Label.FontSize = 14; cb.FontSize = 13;
    end

    % Per-panel legend, top-right corner, symbol + slope value for each line
    legend(ax, [hPoint hYork hOLS hSMA], ...
        {'1\times1^{\circ} grid cell', ...
         sprintf('York: %+.2f', york_b), ...
         sprintf('OLS: %+.2f', ols_b), ...
         sprintf('RMA: %+.2f', sma_b)}, ...
        'Location','northeast', 'FontSize',13, 'Box','on');
end

out_png = fullfile(OUT_DIR,'fig_grid_mld_npp_trend_regression_CBPM_NEMO_v5.png');
exportgraphics(fig, out_png, 'Resolution',200);
fprintf('\nFigure saved: %s\n', out_png);


%% ---------- local functions ----------
function [a,b] = local_weighted_rma(X,Y,W)
% Weighted geometric-mean / reduced-major-axis (SMA) regression.
% Unweighted RMA is b = sign(r)*std(Y)/std(X); this is the
% weighted-moments generalization, using the same cos-lat area weights
% as OLS and (via scaled uncertainties) York, so all three estimators
% share identical spatial weighting.
W = W(:); X = X(:); Y = Y(:);
sw  = sum(W);
mx  = sum(W.*X)/sw;      my  = sum(W.*Y)/sw;
vx  = sum(W.*(X-mx).^2)/sw;
vy  = sum(W.*(Y-my).^2)/sw;
cxy = sum(W.*(X-mx).*(Y-my))/sw;
r = cxy / sqrt(vx*vy);
b = sign(r) * sqrt(vy/vx);
a = my - b*mx;
end
