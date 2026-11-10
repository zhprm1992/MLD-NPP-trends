% Used for: N2 by biome (SI Table S1)
% Original: Manuscript\MLD\manuscript\nature\NCC\R1\code for revision\correlation_bvf_npp\biomes_1M_bvf200_weighted_ot.m (data paths updated to I:\Prj; computation unchanged)
% biomes_1M_bvf200_weighted_ot.m
%
% Stage 2 of 3: biome-split the regridded depth-weighted 0-200m BVF
% (bvf200_<year>.mat, from prj_grid_bvf200_weighted_ot.m), mirroring
% Code\mld_NEMO\biomes_1M_N2_ot.m exactly (fresh inpolygon split against
% the same meanBiom_FAY.shp SEABIOMES shapefile, same idx convention:
% column-major linear index into a (360,180)=[lon,lat] array).
%
% Output: I:\Prj\N2_1M_NEMO\biomes_1M_1dg\bvf200_weighted\<year>\
% bvf200_<biome>.mat -- fields bvf200, lat, lon, idx (matches the
% N2_<biome>.mat naming/field convention of the old pipeline, just with
% bvf200 instead of N2).
%
% Haipeng Zhao
% 2026-07-29

clear
clc

PATH_IN = 'I:\Prj\N2_1M_NEMO\bvf_200_weighted_one';
PATH_OT = 'I:\Prj\N2_1M_NEMO\biomes_1M_1dg\bvf200_weighted';
SEABIOMES = shaperead( ['D:\Users\zhp19\OneDrive - Indiana State University\Research\Projects'...
    '\', 'DINEOF application\biomes\meanBiom_FAY',...
    '\', 'meanBiom_FAY.shp'] );

addpath('D:\Users\zhp19\OneDrive - Indiana State University\Research\Projects\3D bio model\Code\N2_NEMO'); % nameConvert.m

for i = 1 : 20
    yyyy = num2str( i + 2002 );

    d = load( fullfile( PATH_IN, ['bvf200_', yyyy, '.mat'] ) );
    lat = d.lat1; lon = d.lon1;
    [lat_g, lon_g] = meshgrid( lat, lon );
    lon_g = single( lon_g(:) );
    lat_g = single( lat_g(:) );

    for k = 1 : length( SEABIOMES )
        seaBd = SEABIOMES(k);
        idx = inpolygon( lon_g, lat_g, seaBd.X, seaBd.Y );
        idx = find( idx );
        pixelnum = length( idx );
        area = seaBd.Label;
        biomename = nameConvert( area );
        bvf200 = zeros( pixelnum, 12, 'single' );
        for j = 1 : 12
            v1 = d.bvf200G1deg(:, :, j);
            v1 = v1';
            bvf200(:, j) = v1(idx);
        end

        fileName = ['bvf200_', biomename];
        path_ot = fullfile( PATH_OT, yyyy );
        if ~exist( path_ot, 'dir' )
            mkdir( path_ot );
        end
        save( [path_ot, '\', fileName], 'bvf200', 'lat', 'lon', 'idx', '-v7.3' );
    end
    fprintf('%s done\n', yyyy);
end
