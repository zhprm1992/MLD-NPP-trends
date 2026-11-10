% Used for: ensemble NPP by biome
% Original: Code\npp_ensemble\biomes_npp_ensemble_allProducts.m (data paths updated to I:\Prj; computation unchanged)
% Separate ensemble-mean NPP (VGPM, CAFE, CbPM) into biomes, for all 3 MLD
% products (NEMO, ECCO, HYCOM) -- generalizes biomes_npp_ensemble_1M_1dg.m
% (NEMO-only) to loop over products, and writes output using CAFE's
% biome-folder naming (biomes_npp_<product>, no "_ensemble"/"_1M_1dg_mat"
% suffix -- matches npp_cafe_1M_1dg\biomes_npp_<product> exactly).
%
% Haipeng Zhao
% 2026-07-15

clear
clc

products = {'nemo','ecco','hycom'};

SEABIOMES = shaperead( ['D:\Users\zhp19\OneDrive - Indiana State University\Research\Projects'...
    '\', 'DINEOF application\biomes\meanBiom_FAY',...
    '\', 'meanBiom_FAY.shp'] );

PATH_BASE = 'I:\Prj\npp_ensemble_1M_1dg';

for pi = 1:numel(products)
    prod = products{pi};
    fprintf('=== Product: %s ===\n', upper(prod));

    PATH_IN = fullfile(PATH_BASE, ['grid_npp_' prod '_1M_1dg_mat']);
    PATH_OT = fullfile(PATH_BASE, ['biomes_npp_' prod]);

    datafile = dir( fullfile( PATH_IN, '*.mat') );
    for i = 1 : numel( datafile ) % 2003 - 2022
        yyyy = num2str( i + 2002 );

        load( fullfile( datafile(i).folder, datafile(i).name ) ) % nppG1dg, n_models, lat, lon
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
        fprintf('  %s done\n', yyyy);
    end
end

fprintf('\nDone. Output layout mirrors CAFE (npp_cafe_1M_1dg\\biomes_npp_<product>):\n');
for pi = 1:numel(products)
    fprintf('  %s\n', fullfile(PATH_BASE, ['biomes_npp_' products{pi}], '<year>', 'npp_<biome>.mat'));
end
