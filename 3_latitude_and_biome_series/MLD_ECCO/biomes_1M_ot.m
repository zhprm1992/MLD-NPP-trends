% Used for: monthly MLD-ECCO by biome
% Original: Code\mld_ECCO\biomes_1M_ot.m (data paths updated to I:\Prj; computation unchanged)
% Separate global mld into biomes
% Haipeng Zhao
% 9/9/2023
% 9/26/2023
% 10/11/2023 Remove values more than  three standard deviation from mean
% 10/17/2023 Extended to process SST dataset
% 09/25/2024 Applied to POC dataset
% 10/4/2024 Applied to NPP dataset
% 9/5/2025

clear
clc

PATH_IN = 'I:\Prj\MLD_ECCO\one';
PATH_OT = 'I:\Prj\MLD_ECCO\bomes_1M_1dg';
SEABIOMES = shaperead( ['E:\OneDrive - Indiana State University\Research\Projects'...
    '\', 'DINEOF application\biomes\meanBiom_FAY',...
    '\', 'meanBiom_FAY.shp'] );

path = 'I:\Prj\MLD_ARGO\one';
sub_folder1 = [path, '\', '0.03sigmaTheta_ml_', num2str( 2003 ), '.nc'];
files = dir( [sub_folder1, '\', '*.nc'] );
lat = ncread( [files(1).folder, '\', files(1).name], 'LATITUDE' );
lon = ncread( [files(1).folder, '\', files(1).name], 'LONGITUDE' );

datafile = dir( fullfile( PATH_IN, '*.mat') );
for i = 1 : 20
    yyyy = num2str( i + 2002 );

    load( fullfile( datafile(i).folder, datafile(i).name ) )
    [lat_g, lon_g] = meshgrid( lat, lon );
    lon_g = single( lon_g(:) );
    lat_g = single( lat_g(:) );

    for k = 1 : length( SEABIOMES )
        seaBd = SEABIOMES(k);
        idx = inpolygon( lon_g, lat_g, seaBd.X, seaBd.Y ); % find pixels inside of the region   
        idx = find( idx );
        pixelnum = length( idx );
        area = seaBd.Label;
        biomename = nameConvert( area );
        mld = zeros( pixelnum, 12, 'single' );
        for j = 1 : 12
            mld1 = mldG1deg(:, :, j);
            mld1 = mld1';
            mld(:, j) = mld1(idx);            
        end

        fileName = ['mld_', biomename];
        path_ot = fullfile( PATH_OT, yyyy );
        if ~exist( path_ot, 'dir' )
            mkdir( path_ot );
        end
        save( [path_ot, '\', fileName], 'mld', 'lat', 'lon', 'idx', '-v7.3' );
    end
end