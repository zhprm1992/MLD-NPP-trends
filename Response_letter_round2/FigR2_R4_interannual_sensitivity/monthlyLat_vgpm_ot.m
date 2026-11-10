% Used for: monthly latitudinal NPP-VGPM, nppLat_monthly_vgpm.mat (Figs. R2, R4)
% Original: Code\npp_vgpm\monthlyLat_vgpm_ot.m (data paths updated to I:\Prj; computation unchanged)
% Output latitudinal averaged NPP (VGPM) across the global ocean
% Haipeng Zhao
% 7/12/2026

clear
clc

PATH_IN = 'I:\Prj\npp_vgpm_1M_1dg\biomes_npp_vgpm_1M_1dg_mat';
PATH_OT = 'I:\Prj\npp_vgpm_1M_1dg';

d = load( 'm_topo.mat' );
topo = d.topo';

yearfolders = dir( PATH_IN );
yearfolders = yearfolders( [yearfolders.isdir] & ~ismember( {yearfolders.name}, {'.', '..'} ) );
numYears = numel( yearfolders );
startYear = min( str2double( {yearfolders.name} ) );

d = load( fullfile( PATH_IN, num2str( startYear ), 'npp_A EQU.mat' ) );
lon = d.lon;
lat = d.lat;
clear d;

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

mm = repmat( [1 : 12], 1, numYears )';

nppLat = squeeze( nanmean( npp ) );
save( [PATH_OT, '/', 'nppLat_monthly_vgpm'], 'nppLat', 'mm', 'lat' )
