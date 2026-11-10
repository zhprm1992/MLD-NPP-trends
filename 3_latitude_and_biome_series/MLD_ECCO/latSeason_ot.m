% Used for: seasonal latitudinal MLD-ECCO (Fig. 3)
% Original: Code\mld_ECCO\latSeason_ot.m (data paths updated to I:\Prj; computation unchanged)
% Output seasonal latitudinal average
% Haipeng Zhao
% 2/21/2024
% 3/8/2024  
% 1/29/2025
% 7/2/2025 mld
% 7/17/2025 mld

clear
clc

PATH_IN = 'I:\Prj\MLD_ECCO\bomes_1M_1dg';
PATH_OT = 'I:\Prj\MLD_ECCO';
FILENAME_OT = 'mldLat_season';

d = load( 'm_topo.mat' );
topo = d.topo';

d = load( 'I:\Prj\npp_cafe_1M_1dg\biomes_npp_argo\2003\npp_A EQU.mat' );
lon = d.lon;
lat = d.lat;
clear d;

u1 = length( lon ); 
u2 = length( lat );
mld = nan( length( lon ), length( lat ), 20 * 12 );
basin_no = nan( length( lon ), length( lat ) );

datafile = dir( [PATH_IN, '\2003\', '*.mat'] );
filenum = length( datafile );

basin = 0;
for k = 1 : 20
    yyyy = num2str( k + 2002 );
    for j = 1 : filenum 
        regfile = datafile(j);
        filename = regfile.name;
        disp( filename )
    
        d = load( fullfile( PATH_IN, yyyy, filename ) );
        for i = 1 : length( d.idx )
            [I, J] = ind2sub( [u1, u2], d.idx(i) );
            if ( topo(I, J) > -1000 )
                continue
            end
            mld(I, J, (k - 1) * 12 + 1 : k * 12) = d.mld(i, :);
            basin_no(I, J) = basin;
        end
        basin = basin + 1;
    end
    basin = 0;
end

mm = repmat( 1 : 12, 1, 20 )';

idxN = lat > 0;
idxS = lat < 0;

ind1 = mm == 12 | mm == 1 | mm == 2; % Dec, Jan, Feb
mld1 = mld(:, :, ind1);
for i = 1 : 20
    mldwi(:, :, i) = nanmean( mld1(:, :, (i - 1) * 3 + 1 : (i - 1) * 3 + 3), 3);
end

ind2 = mm == 3 | mm == 4 | mm == 5;  % Mar, April, May
mld1 = mld(:, :, ind2);
for i = 1 : 20
    mldsp(:, :, i) = nanmean( mld1(:, :, (i - 1) * 3 + 1 : (i - 1) * 3 + 3), 3);
end

ind3 = mm == 6 | mm == 7 | mm == 8; % June, July, Aug
mld1 = mld(:, :, ind3);
for i = 1 : 20
    mldsu(:, :, i) = nanmean( mld1(:, :, (i - 1) * 3 + 1 : (i - 1) * 3 + 3), 3);
end

ind4 = mm == 9 | mm == 10 | mm == 11; % Sep, Oct, Nov
mld1 = mld(:, :, ind4);
for i = 1 : 20
    mldau(:, :, i) = nanmean( mld1(:, :, (i - 1) * 3 + 1 : (i - 1) * 3 + 3), 3);
end

mldwiS = mldwi(:, idxS, :);
mldwi(:, idxS, :) = mldsu(:, idxS, :);
mldsu(:, idxS, :) = mldwiS;

mldspS = mldsp(:, idxS, :);
mldsp(:, idxS, :) = mldau(:, idxS, :);
mldau(:, idxS, :) = mldspS;

a = squeeze( sum( ~isnan( mldwi ) .* ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
b =  squeeze( sum( ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
c = a ./ b;
mldwiLat = squeeze( nanmean( mldwi ) );
mldwiLat( c < 0.5 ) = nan;

a = squeeze( sum( ~isnan( mldsp ) .* ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
b =  squeeze( sum( ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
c = a ./ b;
mldspLat = squeeze( nanmean( mldsp ) );
mldspLat( c < 0.5 ) = nan;

a = squeeze( sum( ~isnan( mldsu ) .* ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
b =  squeeze( sum( ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
c = a ./ b;
mldsuLat = squeeze( nanmean( mldsu ) );
mldsuLat( c < 0.5 ) = nan;

a = squeeze( sum( ~isnan( mldau ) .* ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
b =  squeeze( sum( ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
c = a ./ b;
mldauLat = squeeze( nanmean( mldau ) );
mldauLat( c < 0.5 ) = nan;

save( [PATH_OT, '/', FILENAME_OT], 'mldwi', 'mldsp', 'mldsu', 'mldau', 'mldwiLat', 'mldspLat', 'mldsuLat', 'mldauLat', 'mm', 'lat' )

