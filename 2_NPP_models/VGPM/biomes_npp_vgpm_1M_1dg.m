% Used for: VGPM NPP split into biomes
% Original: Code\npp_vgpm\biomes_npp_vgpm_1M_1dg.m (data paths updated to I:\Prj; computation unchanged)
% Separate global VGPM NPP into biomes
% Haipeng Zhao
% 7/11/2026

clear
clc

PATH_IN = 'I:\Prj\npp_vgpm_1M_1dg\grid_npp_vgpm_1M_1dg_mat';
PATH_OT = 'I:\Prj\npp_vgpm_1M_1dg\biomes_npp_vgpm_1M_1dg_mat';
SEABIOMES = shaperead( ['D:\Users\zhp19\OneDrive - Indiana State University\Research\Projects'...
    '\', 'DINEOF application\biomes\meanBiom_FAY',...
    '\', 'meanBiom_FAY.shp'] );

datafile = dir( fullfile( PATH_IN, '*.mat') );
for i = 1 : numel( datafile ) % 2003 - 2022
    yyyy = num2str( i + 2002 );

    load( fullfile( datafile(i).folder, datafile(i).name ) ) % nppG1dg, lat, lon
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
        npp = zeros( pixelnum, 12, 'single' );
        for j = 1 : 12
            npp1 = nppG1dg(:, :, j);
            npp1 = npp1';
            npp(:, j) = npp1(idx);
        end

        fileName = ['npp_', biomename];
        path_ot = fullfile( PATH_OT, yyyy );
        if ~exist( path_ot, 'dir' )
            mkdir( path_ot );
        end
        save( [path_ot, '\', fileName], 'npp', 'lat', 'lon', 'idx', '-v7.3' );
    end
end
