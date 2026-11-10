% Used for: index of GLBy0.08 cells in each 1 deg cell (2020-2022, HYCOM grid change)
% Original: Code\mld_HYCOM\proj_GLBy008to1GIdx.m (data paths updated to I:\Prj; computation unchanged)
% Output index of 1 grid from JAMSTEC MLD to which the hycom GLBy008 is projected:
% Haipeng Zhao
% 6/21/2025

clear
clc

%% Target coordinates (ARGO JMSTEC)
fName = 'I:\Prj\MLD_ARGO\one\0.03sigmaTheta_ml_2003.nc\ml_20030101-20030131.nc';
LON = ncread( fName, 'LONGITUDE' );
LAT = ncread( fName, 'LATITUDE' );
[LON, LAT] = meshgrid( LON, LAT );

%%
path_hycom = 'I:\Prj\MLD_1M_HYCOM_GLB\2020\mld_008_1M_2020_01.mat';
load( path_hycom, 'lat' );
load( path_hycom, 'lon' );
lon = mod( lon + 180, 360) - 180;

[lon, lat] = meshgrid( lon, lat );

%% Project variable to 1-deg resolution
gridIdx = nan( length( lon(:) ), 1, 'single' ); % To store which grid the bin belongs to
parfor k = 1 : length( lon(:) )
        gridIdx(k) = search_nearestPixel( lon(k), lat(k), LON(:), LAT(:), 1 );
end

save( 'mldGLBy008to1GIdx', 'gridIdx' )
