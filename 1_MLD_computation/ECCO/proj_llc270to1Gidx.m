% Used for: index of LLC270 cells in each 1 deg cell, mldLLC270to2G_idx.mat
% Original: Code\mld_ECCO\Compute_MLD\proj_llc270to1Gidx.m (data paths updated to I:\Prj; computation unchanged)
% Project LLC270 to 1 deg grid
% Haipeng Zhao
% 8/26/2025

clear
clc

numFaces = 13;
nx = 270;
ny = 270 .* numFaces;
nz = 50;

%% Target coordinates (ARGO JMSTEC)
fName = 'I:\Prj\MLD_ARGO\one\0.03sigmaTheta_ml_2003.nc\ml_20030101-20030131.nc';
LON = ncread( fName, 'LONGITUDE' );
LAT = ncread( fName, 'LATITUDE' );
[LON, LAT] = meshgrid( LON, LAT );

%% Input coordinates (LLC270)
gridDir = 'F:\ECCO-Darwin_extension\grid\';
lat = readbin( [gridDir 'YC.data'], [nx ny],1, 'real*4' );
lon = readbin( [gridDir 'XC.data'], [nx ny],1, 'real*4' );

%% Projection
gridIdx = nan( length( lon(:) ), 1, 'single' ); % To store which grid the bin belongs to
for k = 1 : length( lon(:) )
        gridIdx(k) = search_nearestPixel( lon(k), lat(k), LON(:), LAT(:), 1 );
end

save( 'mldLLC270to2G_idx', 'gridIdx' )