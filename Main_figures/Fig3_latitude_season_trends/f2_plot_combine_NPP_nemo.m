% Used for: Fig. 3a (NEMO row)
% Original: Code\npp_cbpm\f2_plot_combine_NPP_nemo.m (data paths updated to I:\Prj; computation unchanged)
% Plot climatological MLD and MLD trends over the period of 2003 - 2022
% Haipeng Zhao
% 10/16/2025

clear
clc

%% 
lat = ncread( "I:\Prj\ml_20030101-20030131.nc", 'LATITUDE' );
idx = lat <= 60 & lat >= -60;
L = 20; % SET NUMBER OF YEAR

fig = figure('Units', 'inches', 'Position', [1, 1, 15.6, 4]);  % [x, y, width, height]
t = tiledlayout( fig, 1, 8 );
t.TileSpacing = 'tight';
t.Padding = 'compact';
lat = lat(idx);

d1 = load( 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\nppLat_season_nemo.mat' );
d2 = load( 'I:\Prj\MLD_1M_NEMO\mldLat_season.mat' );

cmp = slanCM( 'tab20c' );
cmp = unique( cmp, "rows", 'stable' );
c(1, :) = cmp(7, :);
c(2, :) = cmp(5, :); % Red
c(3, :) = cmp(1, :); % Blue
c(4, :) = cmp(9, :); % Green

c_mld = [31, 119, 180] / 255;   % blue (#1f77b4)
c_npp = [44, 160, 44]  / 255;   % green (#2ca02c)

%% Spring
npp = d1.nppspLat(idx, :);    % Extract range of lat [-60, 60]
mld = d2.mldspLat(idx, :);    % Extract range of lat [-60, 60]

ax1 = axes( t );
h1 = plotClimMLD( ax1, mld, lat,  c_mld );
ax1.XColor = c_mld;

ax2 = axes(t);
h2 = plotClimNPP( ax2, npp, lat, c_npp );
ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
ax2.XColor = c_npp;
ax2.YLabel = [];
ax2.YTick = [];
ax2.Color = 'none';

% legend([h1, h2], {'MLD', 'NPP'}, 'Location', 'southoutside', 'Orientation', 'vertical');
 
nexttile
plotTrendsNPP( npp, mld, lat, 1 )

%% Summer
ax3 = nexttile;
npp = d1.nppsuLat(idx, :);    % Extract range of lat [-60, 60]
mld = d2.mldsuLat(idx, :);    % Extract range of lat [-60, 60]

plotClimMLD( ax3, mld, lat,  c_mld )
ax3.XColor = c_mld;

ax4 = axes(t);
ax4.Layout.Tile=3;
plotClimNPP( ax4, npp, lat, c_npp )
ax4.XAxisLocation = 'top';
ax4.YAxisLocation = 'right';
ax4.XColor = c_npp;
ax4.YLabel = [];
ax4.YTick = [];
ax4.Color = 'none';

nexttile
plotTrendsNPP( npp, mld, lat, 1 )

%%
ax5 = nexttile;
npp = d1.nppauLat(idx, :);    % Extract range of lat [-60, 60]
mld = d2.mldauLat(idx, :);    % Extract range of lat [-60, 60]

plotClimMLD( ax5, mld, lat,  c_mld )
ax5.XColor = c_mld;

ax6 = axes(t);
ax6.Layout.Tile=5;
plotClimNPP( ax6, npp, lat, c_npp )
ax6.XAxisLocation = 'top';
ax6.YAxisLocation = 'right';
ax6.XColor = c_npp;
ax6.YLabel = [];
ax6.YTick = [];
ax6.Color = 'none';

nexttile
plotTrendsNPP( npp, mld, lat, 1 )

%%
ax7 = nexttile;
npp = d1.nppwiLat(idx, :);    % Extract range of lat [-60, 60]
mld = d2.mldwiLat(idx, :);    % Extract range of lat [-60, 60]

plotClimMLD( ax7, mld, lat,  c_mld )
ax7.XColor = c_mld;

ax8 = axes(t);
ax8.Layout.Tile=7;
plotClimNPP( ax8, npp, lat, c_npp )
ax8.XAxisLocation = 'top';
ax8.YAxisLocation = 'right';
ax8.XColor = c_npp;
ax8.YLabel = [];
ax8.YTick = [];
ax8.Color = 'none';

nexttile
plotTrendsNPP( npp, mld, lat, 1 )




