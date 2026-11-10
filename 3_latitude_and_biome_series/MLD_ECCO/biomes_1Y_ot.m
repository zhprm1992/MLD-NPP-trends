% Used for: annual MLD-ECCO by biome (Fig. 2c map)
% Original: Code\mld_ECCO\biomes_1Y_ot.m (data paths updated to I:\Prj; computation unchanged)
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

PATH_IN = 'I:\Prj\MLD_ECCO\one';
FOLDER_OT = 'I:\Prj\MLD_ECCO\biomes_1Y_1dg';
SEABIOMES = shaperead( ['E:\OneDrive - Indiana State University\Research\Projects'...
    '\', 'DINEOF application\biomes\meanBiom_FAY',...
    '\', 'meanBiom_FAY.shp'] );

load( [PATH_IN, '\', 'mld_2003.mat'] )
[lat_g, lon_g] = meshgrid( lat, lon );
lon_g = single( lon_g(:) );
lat_g = single( lat_g(:) );
u = size( mldG1deg );

daynum =  size( mldG1deg, 3 );
for k = 1 : length( SEABIOMES )
    seaBd = SEABIOMES(k);
    idx = inpolygon( lon_g, lat_g, seaBd.X, seaBd.Y ); % find pixels inside of the region   
    idx = find( idx );
    pixelnum = length( idx );
    area = seaBd.Label;
    mld = nan( pixelnum, 20, 'single' );
    for i = 2003 : 2022
        yyyy = num2str( i );
        load( [PATH_IN, '\', 'mld_', yyyy, '.mat'] )
        mld0112 = zeros( pixelnum, daynum, 'single' );
        for j = 1 : daynum
            mld_1M = mldG1deg(:, :, j);
            mld_1M = mld_1M';
            mld0112(:, j) = mld_1M(idx);            
        end
        mld(:, i - 2002) = nanmean( mld0112, 2 );
    end
    fileName = ['mld_', area];
    save( [FOLDER_OT, '\', fileName], 'mld', 'lat', 'lon', 'idx', '-v7.3' );
end




