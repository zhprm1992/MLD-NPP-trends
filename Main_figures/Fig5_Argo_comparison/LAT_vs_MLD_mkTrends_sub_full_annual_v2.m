% Used for: Fig. 5 scatter plots, and the latitudinal bars of Fig. 2
% Original: Code\LAT_vs_MLD_mkTrends_sub_full_annual_v2.m (data paths updated to I:\Prj; computation unchanged)
% Argo mld vs subsampled/fully subsampled model mld 
% Haipeng Zhao
% 4/20/2026

clear
clc

PATH_ARGO = 'I:\Prj\MLD_ARGO\one';
PATH_ECCO = 'I:\Prj\MLD_ECCO\one';
PATH_NEMO = 'I:\Prj\MLD_1M_NEMO\one_d03\mld';
PATH_HYCOM = 'I:\Prj\MLD_HYCOM\one';

color_fullysampled = [0.55 0.25 0.85];
color_fullysampled = [0.7 0.7 0.7];

color_subsampled = [1 0.5 0];
color_subsampled = [0 0 0];

color_nemo= [120, 120, 120] / 255; % NEMO
color_ecco = [230, 159, 0]   / 255; % ECCO
color_hycom = [204, 121, 167]   / 255; % HYCOM

lat = ncread( "I:\Prj\MLD_ARGO\one\0.03sigmaTheta_ml_2003.nc\ml_20030101-20030131.nc", 'LATITUDE' );
idx = lat <= 60 & lat >= -60;
L = 20; % SET NUMBER OF YEAR

for i = 2003 : 2022
    mld_argo = nan( 360, 180, 12, 'single' );
    for j = 1 : 12
            sub_folder1 = [PATH_ARGO, '\', '0.03sigmaTheta_ml_', num2str( i ), '.nc'];
            files = dir( [sub_folder1, '\', '*.nc'] );
            mld_argo(:, :, j) = ncread( [files(j).folder, '\', files(j).name], 'MLD' );
    end
    mld_argo = permute( mld_argo, [2, 1, 3] );

    d_ecco = load(  [PATH_ECCO, '\', 'mld_', num2str( i ), '.mat'] );
    d_nemo = load(  [PATH_NEMO, '\', 'mld_', num2str( i ), '.mat'] );
    d_hycom = load(  [PATH_HYCOM, '\', 'mld_', num2str( i ), '.mat'] );

    mld_nemo_full = d_nemo.mldG1deg;
    mld_ecco_full = d_ecco.mldG1deg;
    mld_hycom_full = d_hycom.mldG1deg;

    mld_nemo_sub = d_nemo.mldG1deg;
    mld_ecco_sub = d_ecco.mldG1deg;
    mld_hycom_sub = d_hycom.mldG1deg;   
    mld_nemo_sub( isnan( mld_argo ) ) = nan;
    mld_ecco_sub( isnan( mld_argo ) ) = nan;
    mld_hycom_sub( isnan( mld_argo ) ) = nan;

    mld_argo_full = nanmean( mld_argo, 3 );
    mld_nemo_full = nanmean( mld_nemo_full, 3 );
    mld_ecco_full = nanmean( mld_ecco_full, 3 );
    mld_hycom_full = nanmean( mld_hycom_full, 3 );

    mld_argo_mothly_sub = nanmean( mld_argo, 3 );
    mld_nemo_sub = nanmean( mld_nemo_sub, 3 );
    mld_ecco_sub = nanmean( mld_ecco_sub, 3 );
    mld_hycom_sub = nanmean( mld_hycom_sub, 3 );

    lat_mld_argo(:, i - 2002) = squeeze( nanmean( mld_argo_full, 2 ) );
    lat_mld_nemo_fulls(:, i - 2002) = squeeze( nanmean( mld_nemo_full, 2 ) );
    lat_mld_ecco_fulls(:, i - 2002) = squeeze( nanmean( mld_ecco_full, 2 ) );
    lat_mld_hycom_fulls(:, i - 2002) = squeeze( nanmean( mld_hycom_full, 2 ) );

    lat_mld_nemo_subs(:, i - 2002) = squeeze( nanmean( mld_nemo_sub, 2 ) );
    lat_mld_ecco_subs(:, i - 2002) = squeeze( nanmean( mld_ecco_sub, 2 ) );
    lat_mld_hycom_subs(:, i - 2002) = squeeze( nanmean( mld_hycom_sub, 2 ) );   
end

% fig = figure;
% t = tiledlayout( fig, 3, 7, 'TileSpacing', 'tight' );
% t.Padding = 'compact';

%% Plot mld_argo vs lat
figure('Position', [100, 100, 150, 350]); % [left, bottom, width, height]

lat_mld_argo = lat_mld_argo(idx, :);    % Extract range of lat [-60, 60]
lat = lat(idx);
p_value = nan( size( lat_mld_argo, 1 ), 1 );
r_square = nan( size( lat_mld_argo, 1 ), 1 );
sen_slope_argo = nan( size( lat_mld_argo, 1 ), 1 );
for i =1 : length( lat )
    y = lat_mld_argo(i, :);
    datain = [(1 : L)', y'];
    [~, ~, ~, sig, ~, ~, ~, sen, ~, plotofslope, ~] = ktaub(datain, 0.9999, 1);       
    if ~isnan(plotofslope)
        sen_slope_argo(i) = single( (plotofslope(2, 2) - plotofslope(1, 2)) / plotofslope(1, 2) * 100 );
        SSres = sum( (y' - plotofslope(:, 2)) .^ 2, 'omitnan' );
        SStot = sum( (y' - nanmean( y' )) .^ 2, 'omitnan' );
        r_square(i) = 1 - SSres / SStot;  % Get R-square
    end
    p_value(i) = sig;
end

hold on
h = barh( lat, sen_slope_argo, 'FaceColor', 'flat' );
c = char( 'S', '', 'N' );
yticklabels(arrayfun(@(x) [num2str(abs(x)), char(176), c(sign(x) + 2)], get( gca, 'YTick' ), 'UniformOutput', false));% set(gca,'YAxisLocation','right')
xlim( [-1.5, 1.5] )
cmp = crameri( 'lajolla', 10 );
rgb = vals2colormap_v2( r_square, cmp, 10, [0, 1] );
colormap( cmp );
colorbar( 'northoutside' )
h.CData = rgb;
set( h, 'Edgecolor', 'flat' );
xlabel( '(% yr^{-1})', 'FontSize', 9 )
% set(gca, 'FontSize', 9 )
xline( -1.5 )
box off
grid off
ax = gca;
xline( ax, ax.XLim(2) )
yline( ax, ax.YLim(2) )

for k = 1 : length( p_value )
    y1 = lat(k);
    x1 = sen_slope_argo(k);
    if p_value(k) < 0.05
        text( 1.0 * sign( x1 ), y1 - 0.5, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
end

% title( 'MLD-Argo', 'FontSize', 5 )

%% Plot mld_nemo vs Lat
figure('Position', [100, 100, 150, 350]); % [left, bottom, width, height]
lat_mld_nemo = lat_mld_nemo_fulls(idx, :);    % Extract range of lat [-60, 60]
p_value = nan( size( lat_mld_nemo, 1 ), 1 );
r_square = nan( size( lat_mld_nemo, 1 ), 1 );
sen_slope = nan( size( lat_mld_nemo, 1 ), 1 );
for i =1 : length( lat )
    y = lat_mld_nemo(i, :);
    datain = [(1 : L)', y'];
    [~, ~, ~, sig, ~, ~, ~, sen, ~, plotofslope, ~] = ktaub(datain, 0.9999, 1);       
    if ~isnan(plotofslope)
        sen_slope(i) = single( (plotofslope(2, 2) - plotofslope(1, 2)) / plotofslope(1, 2) * 100 );
        SSres = sum( (y' - plotofslope(:, 2)) .^ 2, 'omitnan' );
        SStot = sum( (y' - nanmean( y' )) .^ 2, 'omitnan' );
        r_square(i) = 1 - SSres / SStot;  % Get R-square
    end
    p_value(i) = sig;
end

hold on
h = barh( lat, sen_slope, 'FaceColor', 'flat' );
c = char( 'S', '', 'N' );
yticklabels(arrayfun(@(x) [num2str(abs(x)), char(176), c(sign(x) + 2)], get( gca, 'YTick' ), 'UniformOutput', false));% set(gca,'YAxisLocation','right')% set( gca, 'YTick', [] )
xlim( [-1.5, 1.5] )
cmp = crameri( 'lajolla', 10 );
rgb = vals2colormap_v2( r_square, cmp, 10, [0, 1] );
colormap( cmp );
colorbar( 'northoutside' )
h.CData = rgb;
set( h, 'Edgecolor', 'flat' );
xlabel( '(% yr^{-1})', 'FontSize', 9 )

% set(gca, 'FontSize', 9 )
xline( -1.5 )
box off
grid off
ax = gca;
xline( ax, ax.XLim(2) )
yline( ax, ax.YLim(2) )

for k = 1 : length( p_value )
    y1 = lat(k);
    x1 = sen_slope(k);
    if p_value(k) < 0.05
        text( 1.0 * sign( x1 ), y1 - 0.5, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
end
% title( 'MLD-NEMO fully sampled', 'FontSize', 5 )

figure('Position', [100, 100, 150, 350]); % [left, bottom, width, height]
lat_mld_nemo = lat_mld_nemo_subs(idx, :);    % Extract range of lat [-60, 60]
p_value = nan( size( lat_mld_nemo, 1 ), 1 );
r_square = nan( size( lat_mld_nemo, 1 ), 1 );
sen_slope = nan( size( lat_mld_nemo, 1 ), 1 );
for i =1 : length( lat )
    y = lat_mld_nemo(i, :);
    datain = [(1 : L)', y'];
    [~, ~, ~, sig, ~, ~, ~, sen, ~, plotofslope, ~] = ktaub(datain, 0.9999, 1);       
    if ~isnan(plotofslope)
        sen_slope(i) = single( (plotofslope(2, 2) - plotofslope(1, 2)) / plotofslope(1, 2) * 100 );
        SSres = sum( (y' - plotofslope(:, 2)) .^ 2, 'omitnan' );
        SStot = sum( (y' - nanmean( y' )) .^ 2, 'omitnan' );
        r_square(i) = 1 - SSres / SStot;  % Get R-square
    end
    p_value(i) = sig;
end

hold on
h = barh( lat, sen_slope, 'FaceColor', 'flat' );
c = char( 'S', '', 'N' );
yticklabels(arrayfun(@(x) [num2str(abs(x)), char(176), c(sign(x) + 2)], get( gca, 'YTick' ), 'UniformOutput', false));% set(gca,'YAxisLocation','right')% set( gca, 'YTick', [] )
% set(gca,'YAxisLocation','right')
xlim( [-1.5, 1.5] )
cmp = crameri( 'lajolla', 10 );
rgb = vals2colormap_v2( r_square, cmp, 10, [0, 1] );
colormap( cmp );
colorbar( 'northoutside' )
h.CData = rgb;
set( h, 'Edgecolor', 'flat' );
xlabel( '(% yr^{-1})', 'FontSize', 9 )
% set(gca, 'FontSize', 9 )
xline( -1.5 )
box off
grid off
ax = gca;
xline( ax, ax.XLim(2) )
yline( ax, ax.YLim(2) )

for k = 1 : length( p_value )
    y1 = lat(k);
    x1 = sen_slope(k);
    if p_value(k) < 0.05
        text( 1.0 * sign( x1 ), y1 - 0.5, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
end

% title( 'MLD-NEMO subsampled', 'FontSize', 5 )

%% Plot mld_ecco vs Lat
figure('Position', [100, 100, 150, 350]); % [left, bottom, width, height]
lat_mld_ecco = lat_mld_ecco_fulls(idx, :);    % Extract range of lat [-60, 60]
p_value = nan( size( lat_mld_ecco, 1 ), 1 );
r_square = nan( size( lat_mld_ecco, 1 ), 1 );
sen_slope = nan( size( lat_mld_ecco, 1 ), 1 );
for i =1 : length( lat )
    y = lat_mld_ecco(i, :);
    datain = [(1 : L)', y'];
    [~, ~, ~, sig, ~, ~, ~, sen, ~, plotofslope, ~] = ktaub(datain, 0.9999, 1);       
    if ~isnan(plotofslope)
        sen_slope(i) = single( (plotofslope(2, 2) - plotofslope(1, 2)) / plotofslope(1, 2) * 100 );
        SSres = sum( (y' - plotofslope(:, 2)) .^ 2, 'omitnan' );
        SStot = sum( (y' - nanmean( y' )) .^ 2, 'omitnan' );
        r_square(i) = 1 - SSres / SStot;  % Get R-square
    end
    p_value(i) = sig;
end

hold on
h = barh(lat, sen_slope, 'FaceColor', 'flat' );
c = char( 'S', '', 'N' );
yticklabels(arrayfun(@(x) [num2str(abs(x)), char(176), c(sign(x) + 2)], get( gca, 'YTick' ), 'UniformOutput', false));% set(gca,'YAxisLocation','right')
% set(gca,'YAxisLocation','right')
xlim( [-1.5, 1.5] )
cmp = crameri( 'lajolla', 10 );
rgb = vals2colormap_v2( r_square, cmp, 10, [0, 1] );
colormap( cmp );
colorbar( 'northoutside' )
h.CData = rgb;
set( h, 'Edgecolor', 'flat' );
xlabel( '(% yr^{-1})', 'FontSize', 9 )
% set(gca, 'FontSize', 9 )
xline( -1.5 )
box off
grid off
ax = gca;
xline( ax, ax.XLim(2) )
yline( ax, ax.YLim(2) )

for k = 1 : length( p_value )
    y1 = lat(k);
    x1 = sen_slope(k);
    if p_value(k) < 0.05
        text( 1.0 * sign( x1 ), y1 - 0.5, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
end
% title( 'MLD-ECCO fully sampled', 'FontSize', 5 )

figure('Position', [100, 100, 150, 350]); % [left, bottom, width, height]
lat_mld_ecco = lat_mld_ecco_subs(idx, :);    % Extract range of lat [-60, 60]
p_value = nan( size( lat_mld_ecco, 1 ), 1 );
r_square = nan( size( lat_mld_ecco, 1 ), 1 );
sen_slope = nan( size( lat_mld_ecco, 1 ), 1 );
for i =1 : length( lat )
    y = lat_mld_ecco(i, :);
    datain = [(1 : L)', y'];
    [~, ~, ~, sig, ~, ~, ~, sen, ~, plotofslope, ~] = ktaub(datain, 0.9999, 1);       
    if ~isnan(plotofslope)
        sen_slope(i) = single( (plotofslope(2, 2) - plotofslope(1, 2)) / plotofslope(1, 2) * 100 );
        SSres = sum( (y' - plotofslope(:, 2)) .^ 2, 'omitnan' );
        SStot = sum( (y' - nanmean( y' )) .^ 2, 'omitnan' );
        r_square(i) = 1 - SSres / SStot;  % Get R-square
    end
    p_value(i) = sig;
end

hold on
h = barh( lat, sen_slope, 'FaceColor', 'flat' );
c = char( 'S', '', 'N' );
yticklabels(arrayfun(@(x) [num2str(abs(x)), char(176), c(sign(x) + 2)], get( gca, 'YTick' ), 'UniformOutput', false));% set(gca,'YAxisLocation','right')set( gca, 'YTick', [] )
% set(gca,'YAxisLocation','right')
xlim( [-1.5, 1.5] )
cmp = crameri( 'lajolla', 10 );
rgb = vals2colormap_v2( r_square, cmp, 10, [0, 1] );
colormap( cmp );
colorbar( 'northoutside' )
h.CData = rgb;
set( h, 'Edgecolor', 'flat' );
xlabel( '(% yr^{-1})', 'FontSize', 9 )
% set(gca, 'FontSize', 9 )
xline( -1.5 )
box off
grid off
ax = gca;
xline( ax, ax.XLim(2) )
yline( ax, ax.YLim(2) )

for k = 1 : length( p_value )
    y1 = lat(k);
    x1 = sen_slope(k);
    if p_value(k) < 0.05
        text( 1.0 * sign( x1 ), y1 - 0.5, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
end
% title( 'MLD-ECCO subsampled', 'FontSize', 5 )


%% Plot mld_hycom vs Lat
figure('Position', [100, 100, 150, 350]); % [left, bottom, width, height]
lat_mld_hycom = lat_mld_hycom_fulls(idx, :);    % Extract range of lat [-60, 60]
p_value = nan( size( lat_mld_hycom, 1 ), 1 );
r_square = nan( size( lat_mld_hycom, 1 ), 1 );
sen_slope = nan( size( lat_mld_hycom, 1 ), 1 );
for i =1 : length( lat )
    y = lat_mld_hycom(i, :);
    datain = [(1 : L)', y'];
    [~, ~, ~, sig, ~, ~, ~, sen, ~, plotofslope, ~] = ktaub(datain, 0.9999, 1);       
    if ~isnan(plotofslope)
        sen_slope(i) = single( (plotofslope(2, 2) - plotofslope(1, 2)) / plotofslope(1, 2) * 100 );
        SSres = sum( (y' - plotofslope(:, 2)) .^ 2, 'omitnan' );
        SStot = sum( (y' - nanmean( y' )) .^ 2, 'omitnan' );
        r_square(i) = 1 - SSres / SStot;  % Get R-square
    end
    p_value(i) = sig;
end

hold on
h = barh(lat, sen_slope, 'FaceColor', 'flat' );
c = char( 'S', '', 'N' );
yticklabels(arrayfun(@(x) [num2str(abs(x)), char(176), c(sign(x) + 2)], get( gca, 'YTick' ), 'UniformOutput', false));% set(gca,'YAxisLocation','right')

% set(gca,'YAxisLocation','right')
xlim( [-1.5, 1.5] )
cmp = crameri( 'lajolla', 10 );
rgb = vals2colormap_v2( r_square, cmp, 10, [0, 1] );
colormap( cmp );
colorbar( 'northoutside' )
h.CData = rgb;
set( h, 'Edgecolor', 'flat' );
xlabel( '(% yr^{-1})', 'FontSize', 9 )
% set(gca, 'FontSize', 9 )
xline( -1.5 )
box off
grid off
ax = gca;
xline( ax, ax.XLim(2) )
yline( ax, ax.YLim(2) )

for k = 1 : length( p_value )
    y1 = lat(k);
    x1 = sen_slope(k);
    if p_value(k) < 0.05
        text( 1.0 * sign( x1 ), y1 - 0.5, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
end
% title( 'MLD-HYCOM fully sampled ' )
figure('Position', [100, 100, 150, 350]); % [left, bottom, width, height]
lat_mld_hycom = lat_mld_hycom_subs(idx, :);    % Extract range of lat [-60, 60]
p_value = nan( size( lat_mld_hycom, 1 ), 1 );
r_square = nan( size( lat_mld_hycom, 1 ), 1 );
sen_slope = nan( size( lat_mld_hycom, 1 ), 1 );
for i =1 : length( lat )
    y = lat_mld_hycom(i, :);
    datain = [(1 : L)', y'];
    [~, ~, ~, sig, ~, ~, ~, sen, ~, plotofslope, ~] = ktaub(datain, 0.9999, 1);       
    if ~isnan(plotofslope)
        sen_slope(i) = single( (plotofslope(2, 2) - plotofslope(1, 2)) / plotofslope(1, 2) * 100 );
        SSres = sum( (y' - plotofslope(:, 2)) .^ 2, 'omitnan' );
        SStot = sum( (y' - nanmean( y' )) .^ 2, 'omitnan' );
        r_square(i) = 1 - SSres / SStot;  % Get R-square
    end
    p_value(i) = sig;
end

hold on
h = barh( lat, sen_slope, 'FaceColor', 'flat' );
c = char( 'S', '', 'N' );
yticklabels(arrayfun(@(x) [num2str(abs(x)), char(176), c(sign(x) + 2)], get( gca, 'YTick' ), 'UniformOutput', false));% set(gca,'YAxisLocation','right')set( gca, 'YTick', [] )

% set(gca,'YAxisLocation','right')
xlim( [-1.5, 1.5] )
cmp = crameri( 'lajolla', 10 );
rgb = vals2colormap_v2( r_square, cmp, 10, [0, 1] );
colormap( cmp );
colorbar( 'northoutside' )
h.CData = rgb;
set( h, 'Edgecolor', 'flat' );
xlabel( '(% yr^{-1})', 'FontSize', 9 )
% set(gca, 'FontSize', 9 )
xline( -1.5 )
box off
grid off
ax = gca;
xline( ax, ax.XLim(2) )
yline( ax, ax.YLim(2) )

for k = 1 : length( p_value )
    y1 = lat(k);
    x1 = sen_slope(k);
    if p_value(k) < 0.05
        text( 1.0 * sign( x1 ), y1 - 0.5, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
end
% title( 'MLD-HYCOM subsampled ' )

%% Scatter plot NEMO full vs sub sampled
figure('Position', [100, 100, 400, 400]); % [left, bottom, width, height]

lat_mld_nemo_subs = lat_mld_nemo_subs(idx, :);    % Extract range of lat [-60, 60]
p_value_sub = nan( size( lat_mld_nemo_subs, 1 ), 1 );
r_square_sub = nan( size( lat_mld_nemo_subs, 1 ), 1 );
sen_slope_subs = nan( size( lat_mld_nemo_subs, 1 ), 1 );

lat_mld_nemo_fulls = lat_mld_nemo_fulls(idx, :);    % Extract range of lat [-60, 60]
p_value_fulls = nan( size( lat_mld_nemo_fulls, 1 ), 1 );
r_square_fulls = nan( size( lat_mld_nemo_fulls, 1 ), 1 );
sen_slope_fulls = nan( size( lat_mld_nemo_fulls, 1 ), 1 );
for i =1 : length( lat )
    % Subsampling
    y = lat_mld_nemo_subs(i, :);
    datain = [(1 : 20)', y'];
    [~, ~, ~, sig, ~, ~, ~, sen, ~, plotofslope, ~] = ktaub(datain, 0.9999, 1);       
    if ~isnan( plotofslope )
        sen_slope_subs(i) = single( (plotofslope(2, 2) - plotofslope(1, 2)) / plotofslope(1, 2) * 100 );
        SSres = sum( (y' - plotofslope(:, 2)) .^ 2, 'omitnan' );
        SStot = sum( (y' - nanmean( y' )) .^ 2, 'omitnan' );
        r_square_sub(i) = 1 - SSres / SStot;  % Get R-square
    end
    p_value_sub(i) = sig;

    % Full sampling
     y = lat_mld_nemo_fulls(i, :);
    datain = [(1 : 20)', y'];
    [~, ~, ~, sig, ~, ~, ~, sen, ~, plotofslope, ~] = ktaub(datain, 0.9999, 1);       
    if ~isnan( plotofslope )
        sen_slope_fulls(i) = single( (plotofslope(2, 2) - plotofslope(1, 2)) / plotofslope(1, 2) * 100 );
        SSres = sum( (y' - plotofslope(:, 2)) .^ 2, 'omitnan' );
        SStot = sum( (y' - nanmean( y' )) .^ 2, 'omitnan' );
        r_square_fulls(i) = 1 - SSres / SStot;  % Get R-square
    end
    p_value_fulls(i) = sig;   
end

x = sen_slope_argo;  
% ------ Scatter, argo vs sub sample
y1 = sen_slope_subs; 
idx2 = ~isnan( y1 + x ) ; 
y1 = y1(idx2); 
x = x(idx2);

scatter(x, y1, 20, color_subsampled, ...
    'filled', 'MarkerFaceAlpha', 0.7);
hold on

% --- Axis settings ---
xlim([-1.5, 1.5])
ylim([-1.5, 1.5])
caxis([-60 60]);

% xlabel( 'Argo mld trends (% y^{-1})' );
% ylabel( 'NEMO mld trends (% y^{-1})' );

% --- 1:1 reference line ---
plot([-1.5 1.5], [-1.5 1.5], '--', ...
    'LineWidth', 1.0, ...
    'Color', 'k');

% ------Metrics, argo vs sub sample
bias_sub = median(y1 - x);
mae_sub  = median(abs(y1 - x));

mdl = fitlm(x, y1);   % regression: y = slope*x + intercept
slope_sub = mdl.Coefficients.Estimate(2);
intercept = mdl.Coefficients.Estimate(1);
r2_sub = mdl.Rsquared.Ordinary;

%  Regression line 
xx = linspace(-1.5, 1.5, 100);
yy = slope_sub * xx + intercept;

plot(xx, yy, '-', ...
    'LineWidth', 1.5, ...
    'Color', color_subsampled );


% ------- Scatter, argo vs fully sample (orange)
y2 = sen_slope_fulls;
idx2 = ~isnan( y2 + x ) ; 
scatter(x, y2, 20, color_fullysampled - 0.2 );



bias_full = median(y2 - x);
mae_full  = median(abs(y2 - x));

mdl = fitlm(x, y2);   % regression: y = slope*x + intercept
slope_full = mdl.Coefficients.Estimate(2);
intercept = mdl.Coefficients.Estimate(1);
r2_full = mdl.Rsquared.Ordinary;

%  Regression line 
xx = linspace(-1.5, 1.5, 100);
yy = slope_full * xx + intercept;

plot(xx, yy, '-', ...
    'LineWidth', 1.5, ...
    'Color', color_fullysampled );

%------------------------------
% --- Text annotation (top-left) ---
txt = sprintf( ['MAE = %.2f\n', ...
               'Bias = %.2f\n', ...
               'R^2 = %.2f\n', ...
               'Slope = %.2f'], ...
               mae_sub, bias_sub, r2_sub, slope_sub );

text(0.05, 0.95, txt, ...
    'Units', 'normalized', ...
    'VerticalAlignment', 'top', ...
    'FontSize', 15, ...
    'Color', color_subsampled,...
    'BackgroundColor', 'w', ...
    'EdgeColor', 'none');
% --- Styling ---
set(gca, 'FontSize', 12, 'Box', 'on')


%% Scatter plot ECCO full vs sub sampled
figure('Position', [100, 100, 400, 400]); % [left, bottom, width, height]
lat_mld_ecco_subs = lat_mld_ecco_subs(idx, :);    % Extract range of lat [-60, 60]
p_value_sub = nan( size( lat_mld_ecco_subs, 1 ), 1 );
r_square_sub = nan( size( lat_mld_ecco_subs, 1 ), 1 );
sen_slope_subs = nan( size( lat_mld_ecco_subs, 1 ), 1 );

lat_mld_ecco_fulls = lat_mld_ecco_fulls(idx, :);    % Extract range of lat [-60, 60]
p_value_fulls = nan( size( lat_mld_ecco_fulls, 1 ), 1 );
r_square_fulls = nan( size( lat_mld_ecco_fulls, 1 ), 1 );
sen_slope_fulls = nan( size( lat_mld_ecco_fulls, 1 ), 1 );

for i =1 : length( lat )
    % Subsampling
    y = lat_mld_ecco_subs(i, :);
    datain = [(1 : 20)', y'];
    [~, ~, ~, sig, ~, ~, ~, sen, ~, plotofslope, ~] = ktaub(datain, 0.9999, 1);       
    if ~isnan( plotofslope )
        sen_slope_subs(i) = single( (plotofslope(2, 2) - plotofslope(1, 2)) / plotofslope(1, 2) * 100 );
        SSres = sum( (y' - plotofslope(:, 2)) .^ 2, 'omitnan' );
        SStot = sum( (y' - nanmean( y' )) .^ 2, 'omitnan' );
        r_square_sub(i) = 1 - SSres / SStot;  % Get R-square
    end
    p_value_sub(i) = sig;

    % Full sampling
     y = lat_mld_ecco_fulls(i, :);
    datain = [(1 : 20)', y'];
    [~, ~, ~, sig, ~, ~, ~, sen, ~, plotofslope, ~] = ktaub(datain, 0.9999, 1);       
    if ~isnan( plotofslope )
        sen_slope_fulls(i) = single( (plotofslope(2, 2) - plotofslope(1, 2)) / plotofslope(1, 2) * 100 );
        SSres = sum( (y' - plotofslope(:, 2)) .^ 2, 'omitnan' );
        SStot = sum( (y' - nanmean( y' )) .^ 2, 'omitnan' );
        r_square_fulls(i) = 1 - SSres / SStot;  % Get R-square
    end
    p_value_fulls(i) = sig;   
end

x = sen_slope_argo;  
% ------ Scatter, argo vs sub sample
y1 = sen_slope_subs; 
idx2 = ~isnan( y1 + x ) ; 
y1 = y1(idx2); 
x = x(idx2);

scatter(x, y1, 20, color_subsampled, ...
    'filled', 'MarkerFaceAlpha', 0.7);
hold on

% --- Axis settings ---
xlim([-1.5, 1.5])
ylim([-1.5, 1.5])
caxis([-60 60]);

% xlabel( 'Argo mld trends (% y^{-1})' );
% ylabel( 'NEMO mld trends (% y^{-1})' );

% --- 1:1 reference line ---
plot([-1.5 1.5], [-1.5 1.5], '--', ...
    'LineWidth', 1.0, ...
    'Color', 'k');

% ------Metrics, argo vs sub sample
bias_sub = median(y1 - x);
mae_sub  = median(abs(y1 - x));

mdl = fitlm(x, y1);   % regression: y = slope*x + intercept
slope_sub = mdl.Coefficients.Estimate(2);
intercept = mdl.Coefficients.Estimate(1);
r2_sub = mdl.Rsquared.Ordinary;

%  Regression line 
xx = linspace(-1.5, 1.5, 100);
yy = slope_sub * xx + intercept;

plot(xx, yy, '-', ...
    'LineWidth', 1.5, ...
    'Color', color_subsampled );


% ------- Scatter, argo vs fully sample (orange)
y2 = sen_slope_fulls;
idx2 = ~isnan( y2 + x ) ; 
scatter(x, y2, 20, color_fullysampled - 0.2 );



bias_full = median(y2 - x);
mae_full  = median(abs(y2 - x));

mdl = fitlm(x, y2);   % regression: y = slope*x + intercept
slope_full = mdl.Coefficients.Estimate(2);
intercept = mdl.Coefficients.Estimate(1);
r2_full = mdl.Rsquared.Ordinary;

%  Regression line 
xx = linspace(-1.5, 1.5, 100);
yy = slope_full * xx + intercept;

plot(xx, yy, '-', ...
    'LineWidth', 1.5, ...
    'Color', color_fullysampled );

%------------------------------
% --- Text annotation (top-left) ---
txt = sprintf( ['MAE = %.2f\n', ...
               'Bias = %.2f\n', ...
               'R^2 = %.2f\n', ...
               'Slope = %.2f'], ...
               mae_sub, bias_sub, r2_sub, slope_sub );

text(0.05, 0.95, txt, ...
    'Units', 'normalized', ...
    'VerticalAlignment', 'top', ...
    'FontSize', 15, ...
    'Color', color_subsampled,...
    'BackgroundColor', 'w', ...
    'EdgeColor', 'none');

% --- Styling ---
set(gca, 'FontSize', 12, 'Box', 'on')


%% Argo trends vs HYCOM trends
figure('Position', [100, 100, 400, 400]); % [left, bottom, width, height]
lat_mld_hycom_subs = lat_mld_hycom_subs(idx, :);    % Extract range of lat [-60, 60]
p_value_sub = nan( size( lat_mld_hycom_subs, 1 ), 1 );
r_square_sub = nan( size( lat_mld_hycom_subs, 1 ), 1 );
sen_slope_subs = nan( size( lat_mld_hycom_subs, 1 ), 1 );

lat_mld_hycom_fulls = lat_mld_hycom_fulls(idx, :);    % Extract range of lat [-60, 60]
p_value_fulls = nan( size( lat_mld_hycom_fulls, 1 ), 1 );
r_square_fulls = nan( size( lat_mld_hycom_fulls, 1 ), 1 );
sen_slope_fulls = nan( size( lat_mld_hycom_fulls, 1 ), 1 );

for i =1 : length( lat )
    % Subsampling
    y = lat_mld_hycom_subs(i, :);
    datain = [(1 : 20)', y'];
    [~, ~, ~, sig, ~, ~, ~, sen, ~, plotofslope, ~] = ktaub(datain, 0.9999, 1);       
    if ~isnan( plotofslope )
        sen_slope_subs(i) = single( (plotofslope(2, 2) - plotofslope(1, 2)) / plotofslope(1, 2) * 100 );
        SSres = sum( (y' - plotofslope(:, 2)) .^ 2, 'omitnan' );
        SStot = sum( (y' - nanmean( y' )) .^ 2, 'omitnan' );
        r_square_sub(i) = 1 - SSres / SStot;  % Get R-square
    end
    p_value_sub(i) = sig;

    % Full sampling
     y = lat_mld_hycom_fulls(i, :);
    datain = [(1 : 20)', y'];
    [~, ~, ~, sig, ~, ~, ~, sen, ~, plotofslope, ~] = ktaub(datain, 0.9999, 1);       
    if ~isnan( plotofslope )
        sen_slope_fulls(i) = single( (plotofslope(2, 2) - plotofslope(1, 2)) / plotofslope(1, 2) * 100 );
        SSres = sum( (y' - plotofslope(:, 2)) .^ 2, 'omitnan' );
        SStot = sum( (y' - nanmean( y' )) .^ 2, 'omitnan' );
        r_square_fulls(i) = 1 - SSres / SStot;  % Get R-square
    end
    p_value_fulls(i) = sig;   
end

%-------------------------------------------------
x = sen_slope_argo;  
% ------ Scatter, argo vs sub sample
y1 = sen_slope_subs; 
idx2 = ~isnan( y1 + x ) ; 
y1 = y1(idx2); 
x = x(idx2);

scatter(x, y1, 20, color_subsampled, ...
    'filled', 'MarkerFaceAlpha', 0.7);
hold on

% --- Axis settings ---
xlim([-1.5, 1.5])
ylim([-1.5, 1.5])
caxis([-60 60]);

% xlabel( 'Argo mld trends (% y^{-1})' );
% ylabel( 'NEMO mld trends (% y^{-1})' );

% --- 1:1 reference line ---
plot([-1.5 1.5], [-1.5 1.5], '--', ...
    'LineWidth', 1.0, ...
    'Color', 'k');

% ------Metrics, argo vs sub sample
bias_sub = median(y1 - x);
mae_sub  = median(abs(y1 - x));

mdl = fitlm(x, y1);   % regression: y = slope*x + intercept
slope_sub = mdl.Coefficients.Estimate(2);
intercept = mdl.Coefficients.Estimate(1);
r2_sub = mdl.Rsquared.Ordinary;

%  Regression line 
xx = linspace(-1.5, 1.5, 100);
yy = slope_sub * xx + intercept;

plot(xx, yy, '-', ...
    'LineWidth', 1.5, ...
    'Color', color_subsampled );


% ------- Scatter, argo vs fully sample (orange)
y2 = sen_slope_fulls;
idx2 = ~isnan( y2 + x ) ; 
scatter(x, y2, 20, color_fullysampled - 0.2 );


bias_full = median(y2 - x);
mae_full  = median(abs(y2 - x));

mdl = fitlm(x, y2);   % regression: y = slope*x + intercept
slope_full = mdl.Coefficients.Estimate(2);
intercept = mdl.Coefficients.Estimate(1);
r2_full = mdl.Rsquared.Ordinary;

%  Regression line 
xx = linspace(-1.5, 1.5, 100);
yy = slope_full * xx + intercept;

plot(xx, yy, '-', ...
    'LineWidth', 1.5, ...
    'Color', color_fullysampled );

%------------------------------
% --- Text annotation (top-left) ---
txt = sprintf( ['MAE = %.2f\n', ...
               'Bias = %.2f\n', ...
               'R^2 = %.2f\n', ...
               'Slope = %.2f'], ...
               mae_sub, bias_sub, r2_sub, slope_sub );

text(0.05, 0.95, txt, ...
    'Units', 'normalized', ...
    'VerticalAlignment', 'top', ...
    'FontSize', 15, ...
    'Color', color_subsampled,...
    'BackgroundColor', 'w', ...
    'EdgeColor', 'none');

% --- Styling ---
set(gca, 'FontSize', 12, 'Box', 'on' )