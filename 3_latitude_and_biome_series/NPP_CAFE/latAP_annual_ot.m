% Used for: annual latitudinal aph(443), apLat_annual.mat (Extended Data Fig. 1)
% Original: Code\npp_cafe\latAP_annual_ot.m (data paths updated to I:\Prj; computation unchanged)
% Calculate positive and latitudinal mean, input is monthly data splited by
% in biome.
% Haipeng Zhao 
% 10/10/2023
% 3/20/2025
% 5/22/2025
% 6/9/2025

clear
clc 

PATH_IN = 'I:\Prj\npp_cafe_1M_1dg\biomes_npp_nemo';
PATH_OT = 'I:\Prj\npp_cafe_1M_1dg';
FILENAME_OT = 'apLat_annual_nemo';

d = load( [PATH_IN, '\', '2003\npp_A EQU.mat'] );
lon = d.lon;
lat = d.lat;

u1 = length( lon ); 
u2 = length( lat );
ap = nan( length( lon ), length( lat ), 20 );
ap_grid = nan( length( lon ), length( lat ), 20 );

datafile = dir( [PATH_IN, '\2003\', '*.mat'] );
filenum = length( datafile );

for k = 1 : 20
    yyyy = num2str( 2002 + k );
    for j = 1 : filenum 
        regfile = datafile(j);
        filename = regfile.name;
        disp( filename )
    
        d = load( fullfile( PATH_IN, yyyy, filename ) );
        for i = 1 : length( d.idx )
            [I, J] = ind2sub( [u1, u2], d.idx(i) );
            ap_grid(I, J, k) = nanmean( d.ap(i, :), 2 );
        end
    end
end

apLatN = ap_grid(:, 46 : 90, :); % latitude from 1 to 89 degree
apLatS = flip( ap_grid(:, 1 : 45, :), 2 ); % latitude from -1 to -89 degree
apLAT = cat( 1, apLatS, apLatN );
apLatPos = lat(46 : 90);
apLATave = squeeze( mean( apLAT, 1, 'omitnan' ) );

save( [PATH_OT, '/', FILENAME_OT], 'ap_grid', 'apLatN', 'apLatS', 'apLatPos', 'lat', 'lon' )
