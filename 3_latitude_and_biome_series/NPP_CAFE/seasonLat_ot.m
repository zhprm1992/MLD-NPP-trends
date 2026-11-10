% Used for: seasonal latitudinal NPP-CAFE, nppLat_season_<product>.mat (Fig. S2); run once per MLD product (PATH_IN, FILENAME_OT)
% Original: Code\npp_cafe\seasonLat_ot.m (data paths updated to I:\Prj; computation unchanged)
% Output seasonal latitudinal average
% Haipeng Zhao
% 2/21/2024
% 3/8/2024  
% 1/29/2025

clear
clc

PATH_IN = 'I:\Prj\npp_cafe_1M_1dg\biomes_npp_hycom';
PATH_OT = 'I:\Prj\npp_cafe_1M_1dg';
FILENAME_OT = 'nppLat_season_hycom';

d = load( 'm_topo.mat' );
topo = d.topo';

d = load( 'I:\Prj\npp_cafe_1M_1dg\biomes_npp_hycom\2003\npp_A EQU.mat' );
lon = d.lon;
lat = d.lat;
clear d;

u1 = length( lon ); 
u2 = length( lat );
npp = nan( length( lon ), length( lat ), 20 * 12 );
basin_no = nan( length( lon ), length( lat ) );

datafile = dir( [PATH_IN, '\2003\', '*.mat'] );
filenum = length( datafile );

basin = 0;
for k = 1 : 20
    yyyy = num2str( k + 2002 );
    for j = 1 : filenum 
        regfile = datafile(j);
        filename = regfile.name;
%         disp( filename )
        d = load( fullfile( PATH_IN, yyyy, filename ) );
        for i = 1 : length( d.idx )
            [I, J] = ind2sub( [u1, u2], d.idx(i) );
            if ( topo(I, J) > -1000 )
                continue
            end            
            npp(I, J, (k - 1) * 12 + 1 : k * 12) = d.npp(i, :);
            basin_no(I, J) = basin;
        end
        basin = basin + 1;
    end
    basin = 0;
end

mm = repmat([1 : 12], 1, 20)';

ind1 = mm == 12 | mm == 1 | mm == 2; % Dec, Jan, Feb
npp1 = npp(:, :, ind1);
for i = 1 : 20
    nppwi(:, :, i) = nanmean( npp1(:, :, (i - 1) * 3 + 1 : (i - 1) * 3 + 3), 3);
end

ind2 = mm == 3 | mm == 4 | mm == 5;  % Mar, April, May
npp1 = npp(:, :, ind2);
for i = 1 : 20
    nppsp(:, :, i) = nanmean( npp1(:, :, (i - 1) * 3 + 1 : (i - 1) * 3 + 3), 3);
end

ind3 = mm == 6 | mm == 7 | mm == 8; % June, July, Aug
npp1 = npp(:, :, ind3);
for i = 1 : 20
    nppsu(:, :, i) = nanmean( npp1(:, :, (i - 1) * 3 + 1 : (i - 1) * 3 + 3), 3);
end

ind4 = mm == 9 | mm == 10 | mm == 11; % Sep, Oct, Nov
npp1 = npp(:, :, ind4);
for i = 1 : 20
    nppau(:, :, i) = nanmean( npp1(:, :, (i - 1) * 3 + 1 : (i - 1) * 3 + 3), 3);
end

% Switch season in south hemisphere
idxS = lat < 0;

nppwiS = nppwi(:, idxS, :);
nppwi(:, idxS, :) = nppsu(:, idxS, :);
nppsu(:, idxS, :) = nppwiS;

nppspS = nppsp(:, idxS, :);
nppsp(:, idxS, :) = nppau(:, idxS, :);
nppau(:, idxS, :) = nppspS;

a = squeeze( sum( ~isnan( nppwi ) .* ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
b =  squeeze( sum( ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
c = a ./ b;
nppwiLat = squeeze( nanmean( nppwi ) );
nppwiLat( c < 0.5 ) = nan;

a = squeeze( sum( ~isnan( nppsp ) .* ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
b =  squeeze( sum( ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
c = a ./ b;
nppspLat = squeeze( nanmean( nppsp ) );
nppspLat( c < 0.5 ) = nan;


a = squeeze( sum( ~isnan( nppsu ) .* ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
b =  squeeze( sum( ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
c = a ./ b;
nppsuLat = squeeze( nanmean( nppsu ) );
nppsuLat( c < 0.5 ) = nan;

a = squeeze( sum( ~isnan( nppau ) .* ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
b =  squeeze( sum( ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
c = a ./ b;
nppauLat = squeeze( nanmean( nppau ) );
nppauLat( c < 0.5 ) = nan;

save( [PATH_OT, '/', FILENAME_OT], 'nppwi', 'nppsp', 'nppsu', 'nppau', 'nppwiLat', 'nppspLat', 'nppsuLat', 'nppauLat', 'mm', 'lat' )
disp( ['Save successfully: ', FILENAME_OT] )
