% Used for: Extended Data Fig. 1, run once per MLD product (NEMO, ECCO, HYCOM panels)
% Original: Code\npp_cafe\LAT_vs_npp_aph443_group.m (data paths updated to I:\Prj; computation unchanged)
% Plot Lat vs CAFE-NPP and Lat vs aph443

clear
clc

lat = ncread( "I:\Prj\MLD_ARGO\one\0.03sigmaTheta_ml_2003.nc\ml_20030101-20030131.nc", 'LATITUDE' );
idx = lat <= 60 & lat >= -60;
L = 20; % SET NUMBER OF YEAR

fig = figure;
fig.Position = [100, 100, 300, 350];
t = tiledlayout( fig, 2, 2 );
t.Padding = 'compact';
t.TileSpacing = 'compact';
lat = lat(idx);

%% Plot NPP_hycom vs Lat
nexttile( t, [2, 1]) % Full
d = load( 'I:\Prj\npp_cafe_1M_1dg\nppLat_annual_hycom.mat' );
lat_npp = squeeze( nanmean( d.npp_grid ) );
lat_npp = lat_npp(idx, :);    % Extract range of lat [-60, 60]

p_value = nan( size( lat_npp, 1 ), 1 );
r_square = nan( size( lat_npp, 1 ), 1 );
sen_slope_npp = nan( size( lat_npp, 1 ), 1 );
for i =1 : length( lat )
    y = lat_npp(i, :);
    datain = [(1 : L)', y'];
    [~, ~, ~, sig, ~, ~, ~, sen, ~, plotofslope, ~] = ktaub(datain, 0.9999, 1);       
    if ~isnan(plotofslope)
        sen_slope_npp(i) = single( (plotofslope(2, 2) - plotofslope(1, 2)) / plotofslope(1, 2) * 100 );
        SSres = sum( (y' - plotofslope(:, 2)) .^ 2, 'omitnan' );
        SStot = sum( (y' - nanmean( y' )) .^ 2, 'omitnan' );
        r_square(i) = 1 - SSres / SStot;  % Get R-square
    end
    p_value(i) = sig;
end

hold on
h = barh(lat, sen_slope_npp, 'FaceColor', 'flat' );
c = char( 'S', '', 'N' );
yticklabels(arrayfun(@(x) [num2str(abs(x)), char(176), c(sign(x) + 2)], get( gca, 'YTick' ), 'UniformOutput', false));
% set(gca,'YAxisLocation','right')
xlim( [-1.5, 1.5] )
cmp = crameri( 'lajolla', 10 );
% cmp = brewermap( 11, 'Greens' );

rgb = vals2colormap_v2( r_square, cmp, 10, [0, 1] );
colormap( cmp );
% colorbar( 'northoutside' )
h.CData = rgb;
set( h, 'Edgecolor', 'flat' );
xlabel( '(% yr^{-1})', 'FontSize', 9 )
set(gca, 'FontSize', 9 )
xline( -1.5 )
box off
grid off
ax = gca;
xline( ax, ax.XLim(2) )
yline( ax, ax.YLim(2) )

for k = 1 : length( p_value )
    y1 = lat(k);
    x1 = sen_slope_npp(k);
    if p_value(k) < 0.05
        text( 1.0 * sign( x1 ), y1 - 0.5, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
end

title( 'NPP-CAFE' )

%%
ax = nexttile(t, [2, 1]);
d = load( 'I:\Prj\npp_cafe_1M_1dg\apLat_annual.mat' );
lat_ap = squeeze( nanmean( d.ap_grid ) );
lat_ap = lat_ap(idx, :);    % Extract range of lat [-60, 60]

p_value = nan( size( lat_ap, 1 ), 1 );
r_square = nan( size( lat_ap, 1 ), 1 );
sen_slope_a443 = nan( size( lat_ap, 1 ), 1 );

for i =1 : length( lat )
    y = lat_ap(i, :);
    datain = [(1 : L)', y'];
    [~, ~, ~, sig, ~, ~, ~, sen, ~, plotofslope, ~] = ktaub(datain, 0.9999, 1);       
    if ~isnan(plotofslope)
        sen_slope_a443(i) = single( (plotofslope(2, 2) - plotofslope(1, 2)) / plotofslope(1, 2) * 100 );
        SSres = sum( (y' - plotofslope(:, 2)) .^ 2, 'omitnan' );
        SStot = sum( (y' - nanmean( y' )) .^ 2, 'omitnan' );
        r_square(i) = 1 - SSres / SStot;  % Get R-square
    end
    p_value(i) = sig;
end

hold on
h = barh(lat, sen_slope_a443, 'FaceColor', 'flat' );
set( gca, 'YTick', [] )
xlim( [-1.5, 1.5] )
cmp = crameri( 'lajolla', 10 );
rgb = vals2colormap_v2( r_square, cmp, 10, [0, 1] );
colormap( cmp );
colorbar( 'northoutside' )
h.CData = rgb;
set( h, 'Edgecolor', 'flat' );
xlabel( '(% yr^{-1})', 'FontSize', 9 )
set(gca, 'FontSize', 9 )
xline( -1.5 )
box off
grid off
ax = gca;
xline( ax, ax.XLim(2) )
yline( ax, ax.YLim(2) )

for k = 1 : length( p_value )
    y1 = lat(k);
    x1 = sen_slope_a443(k);
    if p_value(k) < 0.05
        text( 1.0 * sign( x1 ), y1 - 0.5, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
end
title( 'aph(443)' )

%%
fig2 = figure;
fig2.Position = [100, 100, 250, 250];
y = sen_slope_a443; x = sen_slope_npp;
idx2 = ~isnan( y + x ) ;
y = y(idx2); x = x(idx2);
% scatter( x, y, 20, 'k', 'filled' )
scatter( x, y, 15, [0.2 0.5 0.8], 'filled', 'MarkerFaceAlpha', 0.8 )
% colormap( cmp );
caxis([0 1]);
xlabel( 'NPP slope (% y-1)' );
ylabel( 'aph(443) slope (% y-1)' );
xlim( [-1.5, 1.5] )
ylim( [-1.5, 1.5] )
hold on
plot( [-1.5 1.5], [-1.5 1.5], '--', 'LineWidth', 1.5, 'Color', [0.6350 0.0780 0.1840] )

% --- Statistics ---
N(1, 1) = length( y );
error_bias(1, 1) = median( (y - x) ./ abs( x ) ) * 100;
c = corrcoef( x, y );
error_r(1, 1) = c(2);
error_mae(1, 1) = median( abs( (y - x) ./ x ) ) * 100;

% Correlation coefficient and p-value
[r_mat, p_mat] = corrcoef( x, y );
r_val = r_mat(1, 2);
p_val = p_mat(1, 2);
r2_val = r_val^2;

% Linear regression slope (using polyfit)
p = polyfit( x, y, 1 );
slope_val = p(1);

% Plot regression line
x_fit = linspace(-1.5, 1.5, 100);
y_fit = polyval(p, x_fit);
plot( x_fit, y_fit, '-', 'LineWidth', 1.5, 'Color', [0.2 0.2 0.2] )

% Format p-value string
if p_val < 0.01
    p_str = 'p < 0.01';
else
    p_str = sprintf('p = %.3f', p_val);
end

% Add text to panel
text_str = sprintf('r = %.2f\nR^2 = %.2f\nSlope = %.2f\n%s', ...
    r_val, r2_val, slope_val, p_str );
text( -1.4, 1.3, text_str, ...
    'FontSize', 10, ...
    'VerticalAlignment', 'top', ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', 'white', ...
    'EdgeColor', 'none' )