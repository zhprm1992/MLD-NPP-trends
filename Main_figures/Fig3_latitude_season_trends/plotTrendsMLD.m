% Used for: Fig. 3 helper function
% Original: Code\npp_cbpm\plotTrendsMLD.m (data paths updated to I:\Prj; computation unchanged)
function h = plotTrendsMLD( lat_mld, lat )
%  plot LAT vs MLD trends
% Haipeng Zhao
% 10/16/2025

% lat_mld = d.mldspLat(idx, :);    % Extract range of lat [-60, 60]

p_value = nan( size( lat_mld, 1 ), 1 );
r_square = nan( size( lat_mld, 1 ), 1 );
sen_slope_mld = nan( size( lat_mld, 1 ), 1 );
for i =1 : length( lat )
    y = lat_mld(i, :);
    datain = [(1 : 20)', y'];
    [~, ~, ~, sig, ~, ~, ~, sen, ~, plotofslope, ~] = ktaub(datain, 0.9999, 1);       
    if ~isnan(plotofslope)
        sen_slope_mld(i) = single( (plotofslope(2, 2) - plotofslope(1, 2)) / abs( plotofslope(1, 2) ) * 100 );
        SSres = sum( (y' - plotofslope(:, 2)) .^ 2, 'omitnan' );
        SStot = sum( (y' - nanmean( y' )) .^ 2, 'omitnan' );
        r_square(i) = 1 - SSres / SStot;  % Get R-square
    end
    p_value(i) = sig;
end

hold on
h = barh(lat, sen_slope_mld, 'FaceColor', 'flat' );
set(gca, 'YLabel', [], 'YTick', []);

% set(gca,'YAxisLocation','right')
xlim( [-2, 2] )
cmp = brewermap( 10, 'Blues' );
rgb = vals2colormap_v2( r_square, cmp, 10, [0, 1] );
colormap( cmp );
colorbar( 'northoutside' )
h.CData = rgb;
set( h, 'Edgecolor', 'flat' );
% xlabel( '(% yr^{-1})', 'FontSize', 9 )
% set(gca, 'FontSize', 9 )
xline( -2 )
box off
grid off
ax = gca;
xline( ax, ax.XLim(2) )
yline( ax, ax.YLim(2) )

for k = 1 : length( p_value )
    y1 = lat(k);
    x1 = sen_slope_mld(k);
    if p_value(k) < 0.05
        text( 1.0 * sign( x1 ), y1 - 0.5, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
end

end