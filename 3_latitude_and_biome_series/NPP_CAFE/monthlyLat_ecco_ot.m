% Used for: monthly latitudinal NPP-CAFE with MLD-ECCO (Fig. S3)
% Original: Code\npp_cafe\monthlyLat_ecco_ot.m (data paths updated to I:\Prj; computation unchanged)
% Output  latitudinal averaged nppa across the global ocean
% Haipeng Zhao 
% 10/15/2023
% 11/09/2023 Extend period 2002 July - 2023
% 1/26/2026

clear
clc 

PATH_IN = 'I:\Prj\npp_cafe_1M_1dg\biomes_npp_ecco';
PATH_OT = 'I:\Prj\npp_cafe_1M_1dg';

d = load( 'm_topo.mat' );
topo = d.topo';

d = load( 'I:\Prj\npp_cafe_1M_1dg\biomes_npp_ecco\2003\npp_A EQU.mat' );
lon = d.lon;
lat = d.lat;

clear d;

u1 = length( lon ); 
u2 = length( lat );
npp = nan( length( lon ), length( lat ), 20 * 12);
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

nppLat = squeeze( nanmean( npp ) );
save( [PATH_OT, '/', 'nppLat_monthly_ecco'], 'nppLat', 'mm', 'lat' )


