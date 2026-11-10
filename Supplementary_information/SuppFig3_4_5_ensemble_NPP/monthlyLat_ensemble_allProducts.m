% Used for: monthly latitudinal ensemble NPP (Supplementary Fig. 3)
% Original: Code\npp_ensemble\monthlyLat_ensemble_allProducts.m (data paths updated to I:\Prj; computation unchanged)
% Output latitudinal averaged NPP (ensemble: VGPM + CAFE + CbPM) across the
% global ocean, for all 3 MLD products (NEMO, ECCO, HYCOM) -- generalizes
% monthlyLat_ensemble_ot.m (NEMO-only) to loop over products, output named
% to match CAFE's convention (nppLat_monthly_<product>.mat, no "_ensemble"
% infix -- the ensemble identity is now implied by folder location, same
% as grid_npp_ensemble_1M_1dg_allProducts.m).
%
% Haipeng Zhao
% 2026-07-15

clear
clc

products = {'nemo','ecco','hycom'};
PATH_BASE = 'I:\Prj\npp_ensemble_1M_1dg';

d = load( 'm_topo.mat' );
topo = d.topo';

for pi = 1:numel(products)
    prod = products{pi};
    fprintf('=== Product: %s ===\n', upper(prod));

    PATH_IN = fullfile(PATH_BASE, ['biomes_npp_' prod]);
    PATH_OT = PATH_BASE;

    yearfolders = dir( PATH_IN );
    yearfolders = yearfolders( [yearfolders.isdir] & ~ismember( {yearfolders.name}, {'.', '..'} ) );
    numYears = numel( yearfolders );
    startYear = min( str2double( {yearfolders.name} ) );

    d2 = load( fullfile( PATH_IN, num2str( startYear ), 'npp_A EQU.mat' ) );
    lon = d2.lon;
    lat = d2.lat;
    clear d2;

    u1 = length( lon );
    u2 = length( lat );
    npp = nan( u1, u2, numYears * 12 );
    basin_no = nan( u1, u2 );

    datafile = dir( [PATH_IN, '\', num2str( startYear ), '\', '*.mat'] );
    filenum = length( datafile );

    basin = 0;
    for k = 1 : numYears
        yyyy = num2str( startYear + k - 1 );
        for j = 1 : filenum
            regfile = datafile(j);
            filename = regfile.name;
            dj = load( fullfile( PATH_IN, yyyy, filename ) );
            for i = 1 : length( dj.idx )
                [I, J] = ind2sub( [u1, u2], dj.idx(i) );
                if ( topo(I, J) > -1000 )
                    continue
                end
                npp(I, J, (k - 1) * 12 + 1 : k * 12) = dj.npp(i, :);
                basin_no(I, J) = basin;
            end
            basin = basin + 1;
        end
        basin = 0;
    end

    mm = repmat( [1 : 12], 1, numYears )';

    nppLat = squeeze( nanmean( npp ) );
    filename_ot = ['nppLat_monthly_' prod];
    save( [PATH_OT, '/', filename_ot], 'nppLat', 'mm', 'lat' )
    fprintf('  saved %s\n', fullfile(PATH_OT, filename_ot));
end
