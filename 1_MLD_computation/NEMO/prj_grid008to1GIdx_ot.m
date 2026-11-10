% Used for: index of 1/12 deg cells in each 1 deg cell of the Argo (JAMSTEC) grid, mldG008to1G_idx.mat
% Original: Code\mld_NEMO\prj_grid008to1GIdx_ot.m (data paths updated to I:\Prj; computation unchanged)
% Project 0.25 to 1 degree
% Haipeng Zhao
% 3/18/2024

clear
clc

%% Target coordinates (ARGO JMSTEC)
fName = 'I:\Prj\MLD_ARGO\one\0.03sigmaTheta_ml_2003.nc\ml_20030101-20030131.nc';
LON = ncread( fName, 'LONGITUDE' );
LAT = ncread( fName, 'LATITUDE' );
[LON, LAT] = meshgrid( LON, LAT );

%% Input coordinates (1/12, NEMO)
lat = single( ncread( 'cmems_mod_glo_phy_my_0.083deg_static.nc', 'latitude' ) );
lon = single( ncread( 'cmems_mod_glo_phy_my_0.083deg_static.nc', 'longitude' ) );
[lon, lat] = meshgrid( lon, lat );

%% Projection
gridIdx = nan( length( lon(:) ), 1, 'single' ); % To store which grid the bin belongs to
parfor k = 1 : length( lon(:) )
        gridIdx(k) = search_nearestPixel( lon(k), lat(k), LON(:), LAT(:), 1 );
end

save( 'mldG008to1G_idx', 'gridIdx' )