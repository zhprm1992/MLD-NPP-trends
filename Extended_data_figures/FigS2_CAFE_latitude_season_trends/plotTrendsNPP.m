% Used for: Fig. S2 helper function
% Original: Code\npp_cafe\plotTrendsNPP.m (data paths updated to I:\Prj; computation unchanged)
function [mld, npp] = plotTrendsNPP( lat_npp, lat_mld, lat, p )
%  Plot LAT vs MLD trends
% Input:
%  - p, if p = 1, plot is on. if p = 0, plot is off.
%
% Output:
% - mld, trends of mld_lat
% - npp, trends of npp_lat
%
% Haipeng Zhao
% 10/16/2025

%% NPP
p_value_npp = nan( size( lat_npp, 1 ), 1 );
r_square_npp = nan( size( lat_npp, 1 ), 1 );
sen_slope_npp = nan( size( lat_npp, 1 ), 1 );
for i =1 : length( lat )
    y = lat_npp(i, :);
    if sum( ~isnan( y ) ) == 0
        continue
    end
    datain = [(1 : 20)', y'];
    [~, ~, ~, sig, ~, ~, ~, sen, ~, plotofslope, ~] = ktaub( datain, 0.9999, 1 );       
    if ~isnan(plotofslope)
        sen_slope_npp(i) = single( (plotofslope(2, 2) - plotofslope(1, 2)) / abs( plotofslope(1, 2) ) * 100 );
        SSres = sum( (y' - plotofslope(:, 2)) .^ 2, 'omitnan' );
        SStot = sum( (y' - nanmean( y' )) .^ 2, 'omitnan' );
        r_square_npp(i) = 1 - SSres / SStot;  % Get R-square
    end
    p_value_npp(i) = sig;
end

cmp_npp = brewermap( 11, 'OrRd' );
cmp_npp = brewermap( 11, 'Greens' );

cmp_npp(1, :) = [];
rgb_npp = vals2colormap_v2( r_square_npp, cmp_npp, 10, [0, 1] );
colormap( cmp_npp );

%% MLD
p_value_mld = nan( size( lat_mld, 1 ), 1 );
r_square_mld = nan( size( lat_mld, 1 ), 1 );
sen_slope_mld = nan( size( lat_mld, 1 ), 1 );
for i =1 : length( lat )
    y = lat_mld(i, :);
    if sum( ~isnan( y ) ) == 0
        continue
    end
    datain = [(1 : 20)', y'];
    [~, ~, ~, sig, ~, ~, ~, sen, ~, plotofslope, ~] = ktaub( datain, 0.9999, 1 );       
    if ~isnan(plotofslope)
        sen_slope_mld(i) = single( (plotofslope(2, 2) - plotofslope(1, 2)) / abs( plotofslope(1, 2) ) * 100 );
        SSres = sum( (y' - plotofslope(:, 2)) .^ 2, 'omitnan' );
        SStot = sum( (y' - nanmean( y' )) .^ 2, 'omitnan' );
        r_square_mld(i) = 1 - SSres / SStot;  % Get R-square
    end
    p_value_mld(i) = sig;
end

cmp_mld = brewermap( 11, 'Blues' );
cmp_mld(1, :) = [];
rgb_mld = vals2colormap_v2( r_square_mld, cmp_mld, 10, [0, 1] );
% colormap( cmp_mld );

%% Plot
if p == 1
    hold on
    for i = 1 : length( lat )
        if sign( sen_slope_mld(i) ) == sign( sen_slope_npp(i) ) % Both on one side
            if  abs( sen_slope_mld(i) ) > abs( sen_slope_npp(i) )
                barh( lat(i), sen_slope_mld(i), 'FaceColor', rgb_mld(i, :), 'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 0.3 );
                barh( lat(i), sen_slope_npp(i), 'FaceColor', rgb_npp(i, :), 'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 0.3 );
            else
                barh( lat(i), sen_slope_npp(i), 'FaceColor', rgb_npp(i, :), 'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 0.3 );
                barh( lat(i), sen_slope_mld(i), 'FaceColor', rgb_mld(i, :), 'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 0.3 );
            end
        else % On opposite side
            barh( lat(i), sen_slope_npp(i), 'FaceColor', rgb_npp(i, :), 'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 0.3 );
            barh( lat(i), sen_slope_mld(i), 'FaceColor', rgb_mld(i, :), 'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 0.3 );
        end
    end
    
    set(gca, 'YLabel', [], 'YTick', []);
    
    % set(gca,'YAxisLocation','right')
    xlim( [-2, 2] )
    % cmp = crameri( 'lajolla', 10 );
    
    colorbar( 'northoutside' )
    xline( -2 )
    box off
    grid off
    ax = gca;
    xline( ax, ax.XLim(2) )
    yline( ax, ax.YLim(2) )
    
    dst1 = max( abs( sen_slope_mld(p_value_mld < 0.05 & abs( sen_slope_mld ) < 1.5) ) );
    dst1 = 1.5;
    dst2 = max( abs( sen_slope_npp(p_value_npp < 0.05 & abs( sen_slope_npp ) < 1.5 ) ) );
    dst2 = 1.5;
    
    for k = 1 : length( p_value_mld )
        y1 = lat(k);
        x1 = sen_slope_mld(k);
        if p_value_mld(k) < 0.05
            text( 1.0 * dst1 * sign( x1 ), y1 - 0.5, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 7, 'Color', [0.1922 0.5098 0.7412] );
        end
    
        x2 = sen_slope_npp(k);
        c_start = [0.9020 0.3333 0.0510];
        c_start = [44, 160, 44]  / 255;

        if p_value_npp(k) < 0.05 && p_value_mld(k) < 0.05 && sign( x2 ) == sign( x1 )
            text( 1.08 * dst1 * sign( x1 ), y1 - 0.5, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 7, 'Color', c_start );
        elseif p_value_npp(k) < 0.05
            text( 1.0 * dst2 * sign( x2 ), y1 - 0.5, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 7, 'Color', c_start );
        end    
    end
end

%%
mld.sen_slope = sen_slope_mld;
mld.rsquare = r_square_mld;
mld.pvalue = p_value_mld;

npp.sen_slope = sen_slope_npp;
npp.rsquare = r_square_npp;
npp.pvalue = p_value_npp;

end