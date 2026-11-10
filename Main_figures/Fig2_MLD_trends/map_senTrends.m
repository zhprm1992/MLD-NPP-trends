% Used for: Fig. 2a-d maps (m per decade), run once per product with PATH = Argo I:\Prj\MLD_ARGO\sen_slopes_yearly, NEMO I:\Prj\MLD_1M_NEMO\sens_slope_yearly\mld, ECCO I:\Prj\MLD_ECCO\sens_slope_yearly, HYCOM I:\Prj\MLD_HYCOM\sen_slopes. Latitudinal bars from Fig5_Argo_comparison\LAT_vs_MLD_mkTrends_sub_full_annual_v2.m (fully sampled panels)
% Original: Code\mld_ARGO_jmstec\map_senTrends.m (data paths updated to I:\Prj; computation unchanged)
% Plot MLD annual trends, pacific centered
% Haipeng Zhao
% 2/16/2023
% 10/12/2023
% 8/16/2025
% 10/22/2025

clear
clc

PATH = 'I:\Prj\MLD_ARGO\sen_slopes_yearly\*.mat';

%% Map coordinates (ARGO JMSTEC)
fName = 'I:\Prj\MLD_ARGO\one\0.03sigmaTheta_ml_2003.nc\ml_20030101-20030131.nc';
lon = ncread( fName, 'LONGITUDE' );
lat = ncread( fName, 'LATITUDE' );

slopeG1deg = nan( length( lon ), length( lat ) );

%%
datafiles = dir( PATH ); 
for i = 1 : length( datafiles )
    d = load( [datafiles(i).folder, '\', datafiles(i).name] );
    slope = nanmean( d.slope, 2 );
    idx = uint32( d.idx );
    for j = 1 : length( idx )
        slopeG1deg(idx(j)) = slope(j);
    end
end

idx = lat < -60 | lat > 60;
slopeG1deg(:, idx) = nan;

figure
caxis( [-20 20] )

% cmp = brewermap(12, '*RdBu');
% cmp = brewermap( 265, '*RdBu');

% cmp( [8, 9], : ) = [0.9546 0.9839 0.8909; 0.9546 0.9839 0.8909];
% cmp( [9, 10], : ) = [0.88 0.88 0.88; 0.88 0.88 0.88];
% colormap( cmp )

cmp = flip( slanCM( 'coolwarm' ) );
colormap( cmp )
% cbarrow

hold on
[lon, idx] = sort( lon );
slopeG1deg = slopeG1deg(idx, :);   % VARIABLE NAMES


%%
% Transform Atlantic-centered to Pacific-centered
lon = [-181; lon];
lon2=mod( lon, 360 );
ind=[212 : 361 1 : 211];
lon2 = lon2(ind);
slopeG1deg = [slopeG1deg(end, :);  slopeG1deg];
slopeG1deg = slopeG1deg(ind, :);
lon2(lon2 < 30)=lon2(lon2 < 30) + 360; %...and subtract 360 to some longitudes

m_proj( 'robinson',   'lat', [-80 80], 'lon', [30, 390]);
m_grid('tickdir', 'in', 'linestyle', 'none', 'fontweight', 'bold', 'fontsize', 18, 'ytick', [-60: 20 : 60]);
m_pcolor( gca, lon2, lat, slopeG1deg' )
m_coast('patch',[0 0 0], 'edgecolor','none');


%% Boundary line
% p1 = m_line( [-180, 180], [0, 0], 'linewi', 2, 'color', 'r', 'linestyle', '--' );
% % p1.Color(4) = 0.50;
% 
% p2 = m_line( [-180, 180], [40, 40], 'linewi', 2, 'color', 'r' );
% % p2.Color(4) = 0.50;
% 
% p3 = m_line( [-180, 180], [-40, -40], 'linewi', 2, 'color', 'r' );
% % p3.Color(4) = 0.50;
% 
% p4 = m_line( [-180, 180], [60, 60], 'linewi', 2, 'color', 'r' );
% p4.Color(4) = 0.50;
% 
% p5 = m_line( [-180, 180], [-60, -60], 'linewi', 2, 'color', 'r' );
% p5.Color(4) = 0.50;