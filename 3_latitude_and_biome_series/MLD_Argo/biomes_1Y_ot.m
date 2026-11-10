% Used for: annual MLD-Argo by biome (Fig. 2a map)
% Original: Code\mld_ARGO_jmstec\biomes_1Y_ot.m (data paths updated to I:\Prj; computation unchanged)
% Separate global NCP into biomes
% Haipeng Zhao
% 9/9/2023
% 9/26/2023
% 10/11/2023 Remove values more than  three standard deviation from mean
% 10/17/2023 Extended to process SST dataset
% 09/25/2024 Applied to POC dataset
% 10/4/2024 Applied to NPP dataset
% 10/15/2024
% 11/24/2024
 % 3/20/2025

clear
clc

PATH_IN = 'I:\Prj\MLD_ARGO\one';
FOLDER_OT = 'I:\Prj\MLD_ARGO\biomes_1Y_1dg';
SEABIOMES = shaperead( ['E:\OneDrive - Indiana State University\Research\Projects'...
    '\', 'DINEOF application\biomes\meanBiom_FAY',...
    '\', 'meanBiom_FAY.shp'] );

fName = 'I:\Prj\MLD_ARGO\one\0.03sigmaTheta_ml_2003.nc\ml_20030101-20030131.nc';
lon = ncread( fName, 'LONGITUDE' );
lat = ncread( fName, 'LATITUDE' );
[lat_g, lon_g] = meshgrid( lat, lon );
lon_g = single( lon_g(:) );
lat_g = single( lat_g(:) );

for k = 1 : length( SEABIOMES )
    seaBd = SEABIOMES(k);
    idx = inpolygon( lon_g, lat_g, seaBd.X, seaBd.Y ); % find pixels inside of the region   
    idx = find( idx );
    pixelnum = length( idx );
    area = seaBd.Label;
    mld = nan( pixelnum, 20, 'single' );
    for i = 2003 : 2022
        yyyy = num2str( i );
        sub_folder1 = [PATH_IN, '\', '0.03sigmaTheta_ml_', yyyy, '.nc'];
        files = dir( [sub_folder1, '\', '*.nc'] );
        mld0112 = zeros( pixelnum, 12, 'single' );
        for j = 1 : 12
            mld_1M = ncread( [files(j).folder, '\', files(j).name], 'MLD' );
            mld0112(:, j) = mld_1M(idx);        
        end
        mld(:, i - 2002) = nanmean( mld0112, 2 );
    end
    fileName = ['mld_', area];
    save( [FOLDER_OT, '\', fileName], 'mld', 'lat', 'lon', 'idx', '-v7.3' );
end




