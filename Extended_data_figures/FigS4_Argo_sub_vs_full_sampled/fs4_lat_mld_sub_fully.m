% Used for: Fig. S4
% Original: Code\fs4_lat_mld_sub_fully.m (data paths updated to I:\Prj; computation unchanged)
% Figure S4
% Haipeng Zhao
% 5/5/2026

clear; clc

%% Settings
PATH_ARGO  = 'I:\Prj\MLD_ARGO\one';
PATH_ECCO  = 'I:\Prj\MLD_ECCO\one';
PATH_NEMO  = 'I:\Prj\MLD_1M_NEMO\one_d03\mld';
PATH_HYCOM = 'I:\Prj\MLD_HYCOM\one';

lat = ncread("I:\Prj\MLD_ARGO\one\0.03sigmaTheta_ml_2003.nc\ml_20030101-20030131.nc", 'LATITUDE');
idx = lat <= 60 & lat >= -60;
L   = 20;

%% Load and compute zonal-mean MLD time series
n_lat = length(lat);
lat_mld_argo        = nan(n_lat, L);
lat_mld_nemo_fulls  = nan(n_lat, L);
lat_mld_ecco_fulls  = nan(n_lat, L);
lat_mld_hycom_fulls = nan(n_lat, L);
lat_mld_nemo_subs   = nan(n_lat, L);
lat_mld_ecco_subs   = nan(n_lat, L);
lat_mld_hycom_subs  = nan(n_lat, L);

for i = 2003:2022
    k = i - 2002;
    mld_argo = nan(360, 180, 12, 'single');
    sub_folder = [PATH_ARGO, '\0.03sigmaTheta_ml_', num2str(i), '.nc'];
    files = dir([sub_folder, '\*.nc']);
    for j = 1:12
        mld_argo(:,:,j) = ncread([files(j).folder, '\', files(j).name], 'MLD');
    end
    mld_argo = permute(mld_argo, [2, 1, 3]);  % -> [lat, lon, month]

    d_ecco  = load([PATH_ECCO,  '\mld_', num2str(i), '.mat']);
    d_nemo  = load([PATH_NEMO,  '\mld_', num2str(i), '.mat']);
    d_hycom = load([PATH_HYCOM, '\mld_', num2str(i), '.mat']);

    nemo  = d_nemo.mldG1deg;
    ecco  = d_ecco.mldG1deg;
    hycom = d_hycom.mldG1deg;

    nemo_sub  = nemo;  nemo_sub( isnan(mld_argo)) = nan;
    ecco_sub  = ecco;  ecco_sub( isnan(mld_argo)) = nan;
    hycom_sub = hycom; hycom_sub(isnan(mld_argo)) = nan;

    zm = @(x) squeeze(mean(mean(x, 3, 'omitnan'), 2, 'omitnan'));
    lat_mld_argo(:,       k) = zm(mld_argo);
    lat_mld_nemo_fulls(:, k) = zm(nemo);
    lat_mld_ecco_fulls(:, k) = zm(ecco);
    lat_mld_hycom_fulls(:,k) = zm(hycom);
    lat_mld_nemo_subs(:,  k) = zm(nemo_sub);
    lat_mld_ecco_subs(:,  k) = zm(ecco_sub);
    lat_mld_hycom_subs(:, k) = zm(hycom_sub);
end

%% Plot 3x3 tiled layout
lat_plt = lat(idx);
fig = figure('Units', 'centimeters', 'Position', [2, 2, 9, 15]);
tl = tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
tl.TileSpacing = 'tight';

panels = {lat_mld_argo(idx,:),        lat_mld_nemo_fulls(idx,:),  lat_mld_nemo_subs(idx,:),  ...
          [],                          lat_mld_ecco_fulls(idx,:),  lat_mld_ecco_subs(idx,:),  ...
          [],                          lat_mld_hycom_fulls(idx,:), lat_mld_hycom_subs(idx,:)};

show_ylbl = logical([1, 0, 0,  0, 1, 0,  0, 1, 0]);  % leftmost non-blank column only
show_xlbl = logical([1, 0, 0,  0, 0, 0,  0, 1, 1]);  % bottom non-blank row per column

for t = 1:9
    nexttile;
    if isempty(panels{t})
        axis off
    else
        plot_sen_slope(lat_plt, panels{t}, L);
        if ~show_ylbl(t), set(gca, 'YTickLabel', []); end
        if ~show_xlbl(t), set(gca, 'XTickLabel', []); xlabel(''); end
    end
end

%% Local function
function plot_sen_slope(lat, data, L)
    n         = length(lat);
    p_value   = nan(n, 1);
    r_square  = nan(n, 1);
    sen_slope = nan(n, 1);

    for i = 1:n
        y = data(i, :);
        [~,~,~,sig,~,~,~,~,~,plotofslope,~] = ktaub([(1:L)', y'], 0.9999, 1);
        if ~isnan(plotofslope)
            sen_slope(i) = single((plotofslope(2,2) - plotofslope(1,2)) / plotofslope(1,2) * 100);
            SSres        = sum((y' - plotofslope(:,2)).^2, 'omitnan');
            SStot        = sum((y' - mean(y, 'omitnan')).^2, 'omitnan');
            r_square(i)  = 1 - SSres / SStot;
        end
        p_value(i) = sig;
    end

    hold on
    h = barh(lat, sen_slope, 'FaceColor', 'flat');
    c = char('S', '', 'N');
    yticklabels(arrayfun(@(x) [num2str(abs(x)), char(176), c(sign(x)+2)], ...
                get(gca,'YTick'), 'UniformOutput', false));
    xlim([-1.5, 1.5])
    cmp = crameri('lajolla', 10);
    colormap(cmp)
    h.CData = vals2colormap_v2(r_square, cmp, 10, [0, 1]);
    set(h, 'Edgecolor', 'flat');
    xlabel('(% yr^{-1})', 'FontSize', 7)
    xline(-1.5); box off; grid off
    ax = gca;
    ax.FontSize = 7;
    xline(ax, ax.XLim(2)); yline(ax, ax.YLim(2))

    for k = find(p_value < 0.05)'
        text(sign(sen_slope(k)), lat(k) - 0.5, '*', ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
end
