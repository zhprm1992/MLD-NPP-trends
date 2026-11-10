% Used for: seasonal latitudinal NPP-CbPM, nppLat_season_<product>.mat (Fig. 3); run once per MLD product (PATH_IN = ...\biomes_npp_<product>, FILENAME_OT = nppLat_season_<product>); only the first block (NPP) is used
% Original: Code\npp_cbpm\seasonLat_ot.m (data paths updated to I:\Prj; computation unchanged)
% Output seasonal latitudinal average
% Haipeng ZhaeccoonppauLat
% 2/21/2024
% 3/8/2024  
% 1/29/2025

clc

%% 1. NPP
clear

PATH_IN = 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\biomes_npp_ecco';
PATH_OT = 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU';
FILENAME_OT = 'nppLat_season_ecco';

d = load( 'm_topo.mat' );
topo = d.topo';

d = load( 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\biomes_npp_ecco\2003\npp_A EQU.mat' );
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

%% 2. MU
clear

PATH_IN = 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\biomes_mu_ecco';
PATH_OT = 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU';
FILENAME_OT = 'muLat_season_ecco';

d = load( 'm_topo.mat' );
topo = d.topo';

d = load( 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\biomes_npp_ecco\2003\npp_A EQU.mat' );
lon = d.lon;
lat = d.lat;
clear d;

u1 = length( lon ); 
u2 = length( lat );
mu = nan( length( lon ), length( lat ), 20 * 12 );
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
            mu(I, J, (k - 1) * 12 + 1 : k * 12) = d.mu(i, :);
            basin_no(I, J) = basin;
        end
        basin = basin + 1;
    end
    basin = 0;
end

mm = repmat([1 : 12], 1, 20)';

ind1 = mm == 12 | mm == 1 | mm == 2; % Dec, Jan, Feb
mu1 = mu(:, :, ind1);
for i = 1 : 20
    muwi(:, :, i) = nanmean( mu1(:, :, (i - 1) * 3 + 1 : (i - 1) * 3 + 3), 3);
end

ind2 = mm == 3 | mm == 4 | mm == 5;  % Mar, April, May
mu1 = mu(:, :, ind2);
for i = 1 : 20
    musp(:, :, i) = nanmean( mu1(:, :, (i - 1) * 3 + 1 : (i - 1) * 3 + 3), 3);
end

ind3 = mm == 6 | mm == 7 | mm == 8; % June, July, Aug
mu1 = mu(:, :, ind3);
for i = 1 : 20
    musu(:, :, i) = nanmean( mu1(:, :, (i - 1) * 3 + 1 : (i - 1) * 3 + 3), 3);
end

ind4 = mm == 9 | mm == 10 | mm == 11; % Sep, Oct, Nov
mu1 = mu(:, :, ind4);
for i = 1 : 20
    muau(:, :, i) = nanmean( mu1(:, :, (i - 1) * 3 + 1 : (i - 1) * 3 + 3), 3);
end

idxN = lat > 0;
idxS = lat < 0;

muwiS = muwi(:, idxS, :);
muwi(:, idxS, :) = musu(:, idxS, :);
musu(:, idxS, :) = muwiS;

muspS = musp(:, idxS, :);
musp(:, idxS, :) = muau(:, idxS, :);
muau(:, idxS, :) = muspS;

a = squeeze( sum( ~isnan( muwi ) .* ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
b =  squeeze( sum( ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
c = a ./ b;
muwiLat = squeeze( nanmean( muwi ) );
muwiLat( c < 0.5 ) = nan;

a = squeeze( sum( ~isnan( musp ) .* ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
b =  squeeze( sum( ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
c = a ./ b;
muspLat = squeeze( nanmean( musp ) );
muspLat( c < 0.5 ) = nan;

a = squeeze( sum( ~isnan( musu ) .* ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
b =  squeeze( sum( ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
c = a ./ b;
musuLat = squeeze( nanmean( musu ) );
musuLat( c < 0.5 ) = nan;

a = squeeze( sum( ~isnan( muau ) .* ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
b =  squeeze( sum( ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
c = a ./ b;
muauLat = squeeze( nanmean( muau ) );
muauLat( c < 0.5 ) = nan;

save( [PATH_OT, '/', FILENAME_OT], 'muwi', 'musp', 'musu', 'muau', 'muwiLat', 'muspLat', 'musuLat', 'muauLat', 'mm', 'lat' )
disp( ['Save successfully: ', FILENAME_OT] )

%% 3. CARBON
clear

PATH_IN = 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\biomes_carbon_ecco';
PATH_OT = 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU';
FILENAME_OT = 'carbonLat_season_ecco';

d = load( 'm_topo.mat' );
topo = d.topo';

d = load( 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\biomes_npp_ecco\2003\npp_A EQU.mat' );
lon = d.lon;
lat = d.lat;
clear d;

u1 = length( lon ); 
u2 = length( lat );
carbon = nan( length( lon ), length( lat ), 20 * 12 );
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
            carbon(I, J, (k - 1) * 12 + 1 : k * 12) = d.carbon(i, :);
            basin_no(I, J) = basin;
        end
        basin = basin + 1;
    end
    basin = 0;
end

mm = repmat([1 : 12], 1, 20)';

ind1 = mm == 12 | mm == 1 | mm == 2; % Dec, Jan, Feb
carbon1 = carbon(:, :, ind1);
for i = 1 : 20
    carbonwi(:, :, i) = nanmean( carbon1(:, :, (i - 1) * 3 + 1 : (i - 1) * 3 + 3), 3);
end

ind2 = mm == 3 | mm == 4 | mm == 5;  % Mar, April, May
carbon1 = carbon(:, :, ind2);
for i = 1 : 20
    carbonsp(:, :, i) = nanmean( carbon1(:, :, (i - 1) * 3 + 1 : (i - 1) * 3 + 3), 3);
end

ind3 = mm == 6 | mm == 7 | mm == 8; % June, July, Aug
carbon1 = carbon(:, :, ind3);
for i = 1 : 20
    carbonsu(:, :, i) = nanmean( carbon1(:, :, (i - 1) * 3 + 1 : (i - 1) * 3 + 3), 3);
end

ind4 = mm == 9 | mm == 10 | mm == 11; % Sep, Oct, Nov
carbon1 = carbon(:, :, ind4);
for i = 1 : 20
    carbonau(:, :, i) = nanmean( carbon1(:, :, (i - 1) * 3 + 1 : (i - 1) * 3 + 3), 3);
end

idxN = lat > 0;
idxS = lat < 0;

carbonwiS = carbonwi(:, idxS, :);
carbonwi(:, idxS, :) = carbonsu(:, idxS, :);
carbonsu(:, idxS, :) = carbonwiS;

carbonspS = carbonsp(:, idxS, :);
carbonsp(:, idxS, :) = carbonau(:, idxS, :);
carbonau(:, idxS, :) = carbonspS;

a = squeeze( sum( ~isnan( carbonwi ) .* ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
b =  squeeze( sum( ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
c = a ./ b;
carbonwiLat = squeeze( nanmean( carbonwi ) );
carbonwiLat( c < 0.5 ) = nan;

a = squeeze( sum( ~isnan( carbonsp ) .* ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
b =  squeeze( sum( ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
c = a ./ b;
carbonspLat = squeeze( nanmean( carbonsp ) );
carbonspLat( c < 0.5 ) = nan;

a = squeeze( sum( ~isnan( carbonsu ) .* ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
b =  squeeze( sum( ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
c = a ./ b;
carbonsuLat = squeeze( nanmean( carbonsu ) );
carbonsuLat( c < 0.5 ) = nan;

a = squeeze( sum( ~isnan( carbonau ) .* ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
b =  squeeze( sum( ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
c = a ./ b;
carbonauLat = squeeze( nanmean( carbonau ) );
carbonauLat( c < 0.5 ) = nan;

save( [PATH_OT, '/', FILENAME_OT], 'carbonwi', 'carbonsp', 'carbonsu', 'carbonau', 'carbonwiLat', 'carbonspLat', 'carbonsuLat', 'carbonauLat', 'mm', 'lat' )
disp( ['Save successfully: ', FILENAME_OT] )



%% 4. IgFUN
clear

PATH_IN = 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\biomes_IgFun_ecco';
PATH_OT = 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU';
FILENAME_OT = 'IgFunLat_season_ecco';

d = load( 'm_topo.mat' );
topo = d.topo';

d = load( 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\biomes_npp_ecco\2003\npp_A EQU.mat' );
lon = d.lon;
lat = d.lat;
clear d;

u1 = length( lon ); 
u2 = length( lat );
IgFun = nan( length( lon ), length( lat ), 20 * 12 );
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
            IgFun(I, J, (k - 1) * 12 + 1 : k * 12) = d.IgFun(i, :);
            basin_no(I, J) = basin;
        end
        basin = basin + 1;
    end
    basin = 0;
end

mm = repmat([1 : 12], 1, 20)';

ind1 = mm == 12 | mm == 1 | mm == 2; % Dec, Jan, Feb
IgFun1 = IgFun(:, :, ind1);
for i = 1 : 20
    IgFunwi(:, :, i) = nanmean( IgFun1(:, :, (i - 1) * 3 + 1 : (i - 1) * 3 + 3), 3);
end

ind2 = mm == 3 | mm == 4 | mm == 5;  % Mar, April, May
IgFun1 = IgFun(:, :, ind2);
for i = 1 : 20
    IgFunsp(:, :, i) = nanmean( IgFun1(:, :, (i - 1) * 3 + 1 : (i - 1) * 3 + 3), 3);
end

ind3 = mm == 6 | mm == 7 | mm == 8; % June, July, Aug
IgFun1 = IgFun(:, :, ind3);
for i = 1 : 20
    IgFunsu(:, :, i) = nanmean( IgFun1(:, :, (i - 1) * 3 + 1 : (i - 1) * 3 + 3), 3);
end

ind4 = mm == 9 | mm == 10 | mm == 11; % Sep, Oct, Nov
IgFun1 = IgFun(:, :, ind4);
for i = 1 : 20
    IgFunau(:, :, i) = nanmean( IgFun1(:, :, (i - 1) * 3 + 1 : (i - 1) * 3 + 3), 3);
end

idxN = lat > 0;
idxS = lat < 0;

IgFunwiS = IgFunwi(:, idxS, :);
IgFunwi(:, idxS, :) = IgFunsu(:, idxS, :);
IgFunsu(:, idxS, :) = IgFunwiS;

IgFunspS = IgFunsp(:, idxS, :);
IgFunsp(:, idxS, :) = IgFunau(:, idxS, :);
IgFunau(:, idxS, :) = IgFunspS;

a = squeeze( sum( ~isnan( IgFunwi ) .* ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
b =  squeeze( sum( ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
c = a ./ b;
IgFunwiLat = squeeze( nanmean( IgFunwi ) );
IgFunwiLat( c < 0.5 ) = nan;

a = squeeze( sum( ~isnan( IgFunsp ) .* ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
b =  squeeze( sum( ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
c = a ./ b;
IgFunspLat = squeeze( nanmean( IgFunsp ) );
IgFunspLat( c < 0.5 ) = nan;

a = squeeze( sum( ~isnan( IgFunsu ) .* ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
b =  squeeze( sum( ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
c = a ./ b;
IgFunsuLat = squeeze( nanmean( IgFunsu ) );
IgFunsuLat( c < 0.5 ) = nan;

a = squeeze( sum( ~isnan( IgFunau ) .* ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
b =  squeeze( sum( ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
c = a ./ b;
IgFunauLat = squeeze( nanmean( IgFunau ) );
IgFunauLat( c < 0.5 ) = nan;

save( [PATH_OT, '/', FILENAME_OT], 'IgFunwi', 'IgFunsp', 'IgFunsu', 'IgFunau', 'IgFunwiLat', 'IgFunspLat', 'IgFunsuLat', 'IgFunauLat', 'mm', 'lat' )
disp( ['Save successfully: ', FILENAME_OT] )


%% 5. ntFUN
clear

PATH_IN = 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\biomes_ntFun_ecco';
PATH_OT = 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU';
FILENAME_OT = 'ntFunLat_season_ecco';

d = load( 'm_topo.mat' );
topo = d.topo';

d = load( 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\biomes_npp_ecco\2003\npp_A EQU.mat' );
lon = d.lon;
lat = d.lat;
clear d;

u1 = length( lon ); 
u2 = length( lat );
ntFun = nan( length( lon ), length( lat ), 20 * 12 );
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
            ntFun(I, J, (k - 1) * 12 + 1 : k * 12) = d.ntFun(i, :);
            basin_no(I, J) = basin;
        end
        basin = basin + 1;
    end
    basin = 0;
end

mm = repmat([1 : 12], 1, 20)';

ind1 = mm == 12 | mm == 1 | mm == 2; % Dec, Jan, Feb
ntFun1 = ntFun(:, :, ind1);
for i = 1 : 20
    ntFunwi(:, :, i) = nanmean( ntFun1(:, :, (i - 1) * 3 + 1 : (i - 1) * 3 + 3), 3);
end

ind2 = mm == 3 | mm == 4 | mm == 5;  % Mar, April, May
ntFun1 = ntFun(:, :, ind2);
for i = 1 : 20
    ntFunsp(:, :, i) = nanmean( ntFun1(:, :, (i - 1) * 3 + 1 : (i - 1) * 3 + 3), 3);
end

ind3 = mm == 6 | mm == 7 | mm == 8; % June, July, Aug
ntFun1 = ntFun(:, :, ind3);
for i = 1 : 20
    ntFunsu(:, :, i) = nanmean( ntFun1(:, :, (i - 1) * 3 + 1 : (i - 1) * 3 + 3), 3);
end

ind4 = mm == 9 | mm == 10 | mm == 11; % Sep, Oct, Nov
ntFun1 = ntFun(:, :, ind4);
for i = 1 : 20
    ntFunau(:, :, i) = nanmean( ntFun1(:, :, (i - 1) * 3 + 1 : (i - 1) * 3 + 3), 3);
end

idxN = lat > 0;
idxS = lat < 0;

ntFunwiS = ntFunwi(:, idxS, :);
ntFunwi(:, idxS, :) = ntFunsu(:, idxS, :);
ntFunsu(:, idxS, :) = ntFunwiS;

ntFunspS = ntFunsp(:, idxS, :);
ntFunsp(:, idxS, :) = ntFunau(:, idxS, :);
ntFunau(:, idxS, :) = ntFunspS;

a = squeeze( sum( ~isnan( ntFunwi ) .* ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
b =  squeeze( sum( ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
c = a ./ b;
ntFunwiLat = squeeze( nanmean( ntFunwi ) );
ntFunwiLat( c < 0.5 ) = nan;

a = squeeze( sum( ~isnan( ntFunsp ) .* ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
b =  squeeze( sum( ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
c = a ./ b;
ntFunspLat = squeeze( nanmean( ntFunsp ) );
ntFunspLat( c < 0.5 ) = nan;

a = squeeze( sum( ~isnan( ntFunsu ) .* ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
b =  squeeze( sum( ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
c = a ./ b;
ntFunsuLat = squeeze( nanmean( ntFunsu ) );
ntFunsuLat( c < 0.5 ) = nan;

a = squeeze( sum( ~isnan( ntFunau ) .* ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
b =  squeeze( sum( ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
c = a ./ b;
ntFunauLat = squeeze( nanmean( ntFunau ) );
ntFunauLat( c < 0.5 ) = nan;

save( [PATH_OT, '/', FILENAME_OT], 'ntFunwi', 'ntFunsp', 'ntFunsu', 'ntFunau', 'ntFunwiLat', 'ntFunspLat', 'ntFunsuLat', 'ntFunauLat', 'mm', 'lat' )
disp( ['Save successfully: ', FILENAME_OT] )



%% 6. CHL/CARBON
clear

PATH_IN = 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\biomes_chlCarbon_ecco';
PATH_OT = 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU';
FILENAME_OT = 'chlCarbonLat_season_ecco';

d = load( 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\biomes_npp_ecco\2003\npp_A EQU.mat' );
lon = d.lon;
lat = d.lat;
clear d;

d = load( 'm_topo.mat' );
topo = d.topo';

u1 = length( lon ); 
u2 = length( lat );
chlCarbon = nan( length( lon ), length( lat ), 20 * 12 );
basin_no = nan( length( lon ), length( lat ) );

datafile = dir( [PATH_IN, '\2003\', '*.mat'] );
filenum = length( datafile );

basin = 0;
for k = 1 : 20
    yyyy = num2str( k + 2002 );
    for j = 1 : filenum 
        regfile = datafile(j);
        filename = regfile.name;
        d = load( fullfile( PATH_IN, yyyy, filename ) );
        for i = 1 : length( d.idx )
            [I, J] = ind2sub( [u1, u2], d.idx(i) );
            if ( topo(I, J) > -1000 )
                continue
            end            
            chlCarbon(I, J, (k - 1) * 12 + 1 : k * 12) = d.chlCarbon(i, :);
            basin_no(I, J) = basin;
        end
        basin = basin + 1;
    end
    basin = 0;
end

mm = repmat([1 : 12], 1, 20)';

ind1 = mm == 12 | mm == 1 | mm == 2; % Dec, Jan, Feb
chlCarbon1 = chlCarbon(:, :, ind1);
for i = 1 : 20
    chlCarbonwi(:, :, i) = nanmean( chlCarbon1(:, :, (i - 1) * 3 + 1 : (i - 1) * 3 + 3), 3);
end

ind2 = mm == 3 | mm == 4 | mm == 5;  % Mar, April, May
chlCarbon1 = chlCarbon(:, :, ind2);
for i = 1 : 20
    chlCarbonsp(:, :, i) = nanmean( chlCarbon1(:, :, (i - 1) * 3 + 1 : (i - 1) * 3 + 3), 3);
end

ind3 = mm == 6 | mm == 7 | mm == 8; % June, July, Aug
chlCarbon1 = chlCarbon(:, :, ind3);
for i = 1 : 20
    chlCarbonsu(:, :, i) = nanmean( chlCarbon1(:, :, (i - 1) * 3 + 1 : (i - 1) * 3 + 3), 3);
end

ind4 = mm == 9 | mm == 10 | mm == 11; % Sep, Oct, Nov
chlCarbon1 = chlCarbon(:, :, ind4);
for i = 1 : 20
    chlCarbonau(:, :, i) = nanmean( chlCarbon1(:, :, (i - 1) * 3 + 1 : (i - 1) * 3 + 3), 3);
end

idxN = lat > 0;
idxS = lat < 0;

chlCarbonwiS = chlCarbonwi(:, idxS, :);
chlCarbonwi(:, idxS, :) = chlCarbonsu(:, idxS, :);
chlCarbonsu(:, idxS, :) = chlCarbonwiS;

chlCarbonspS = chlCarbonsp(:, idxS, :);
chlCarbonsp(:, idxS, :) = chlCarbonau(:, idxS, :);
chlCarbonau(:, idxS, :) = chlCarbonspS;

a = squeeze( sum( ~isnan( chlCarbonwi ) .* ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
b =  squeeze( sum( ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
c = a ./ b;
chlCarbonwiLat = squeeze( nanmean( chlCarbonwi ) );
chlCarbonwiLat( c < 0.5 ) = nan;

a = squeeze( sum( ~isnan( chlCarbonsp ) .* ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
b =  squeeze( sum( ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
c = a ./ b;
chlCarbonspLat = squeeze( nanmean( chlCarbonsp ) );
chlCarbonspLat( c < 0.5 ) = nan;

a = squeeze( sum( ~isnan( chlCarbonsu ) .* ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
b =  squeeze( sum( ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
c = a ./ b;
chlCarbonsuLat = squeeze( nanmean( chlCarbonsu ) );
chlCarbonsuLat( c < 0.5 ) = nan;

a = squeeze( sum( ~isnan( chlCarbonau ) .* ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
b =  squeeze( sum( ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
c = a ./ b;
chlCarbonauLat = squeeze( nanmean( chlCarbonau ) );
chlCarbonauLat( c < 0.5 ) = nan;

save( [PATH_OT, '/', FILENAME_OT], 'chlCarbonwi', 'chlCarbonsp', 'chlCarbonsu', 'chlCarbonau', 'chlCarbonwiLat', 'chlCarbonspLat', 'chlCarbonsuLat', 'chlCarbonauLat', 'mm', 'lat' )
disp( ['Save successfully: ', FILENAME_OT] )


%% 7. CHL/CARBON max

clear

PATH_IN = 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\biomes_chlCarbonMx_ecco';
PATH_OT = 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU';
FILENAME_OT = 'chlCarbonMxLat_season_ecco';

d = load( 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\biomes_npp_ecco\2003\npp_A EQU.mat' );
lon = d.lon;
lat = d.lat;
clear d;

d = load( 'm_topo.mat' );
topo = d.topo';

u1 = length( lon ); 
u2 = length( lat );
chlCarbonMx = nan( length( lon ), length( lat ), 20 * 12 );
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
            chlCarbonMx(I, J, (k - 1) * 12 + 1 : k * 12) = d.chlCarbonMx(i, :);
            basin_no(I, J) = basin;
        end
        basin = basin + 1;
    end
    basin = 0;
end

mm = repmat([1 : 12], 1, 20)';

ind1 = mm == 12 | mm == 1 | mm == 2; % Dec, Jan, Feb
chlCarbonMx1 = chlCarbonMx(:, :, ind1);
for i = 1 : 20
    chlCarbonMxwi(:, :, i) = nanmean( chlCarbonMx1(:, :, (i - 1) * 3 + 1 : (i - 1) * 3 + 3), 3);
end

ind2 = mm == 3 | mm == 4 | mm == 5;  % Mar, April, May
chlCarbonMx1 = chlCarbonMx(:, :, ind2);
for i = 1 : 20
    chlCarbonMxsp(:, :, i) = nanmean( chlCarbonMx1(:, :, (i - 1) * 3 + 1 : (i - 1) * 3 + 3), 3);
end

ind3 = mm == 6 | mm == 7 | mm == 8; % June, July, Aug
chlCarbonMx1 = chlCarbonMx(:, :, ind3);
for i = 1 : 20
    chlCarbonMxsu(:, :, i) = nanmean( chlCarbonMx1(:, :, (i - 1) * 3 + 1 : (i - 1) * 3 + 3), 3);
end

ind4 = mm == 9 | mm == 10 | mm == 11; % Sep, Oct, Nov
chlCarbonMx1 = chlCarbonMx(:, :, ind4);
for i = 1 : 20
    chlCarbonMxau(:, :, i) = nanmean( chlCarbonMx1(:, :, (i - 1) * 3 + 1 : (i - 1) * 3 + 3), 3);
end

idxN = lat > 0;
idxS = lat < 0;

chlCarbonMxwiS = chlCarbonMxwi(:, idxS, :);
chlCarbonMxwi(:, idxS, :) = chlCarbonMxsu(:, idxS, :);
chlCarbonMxsu(:, idxS, :) = chlCarbonMxwiS;

chlCarbonMxspS = chlCarbonMxsp(:, idxS, :);
chlCarbonMxsp(:, idxS, :) = chlCarbonMxau(:, idxS, :);
chlCarbonMxau(:, idxS, :) = chlCarbonMxspS;

a = squeeze( sum( ~isnan( chlCarbonMxwi ) .* ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
b =  squeeze( sum( ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
c = a ./ b;
chlCarbonMxwiLat = squeeze( nanmean( chlCarbonMxwi ) );
chlCarbonMxwiLat( c < 0.5 ) = nan;

a = squeeze( sum( ~isnan( chlCarbonMxsp ) .* ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
b =  squeeze( sum( ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
c = a ./ b;
chlCarbonMxspLat = squeeze( nanmean( chlCarbonMxsp ) );
chlCarbonMxspLat( c < 0.5 ) = nan;

a = squeeze( sum( ~isnan( chlCarbonMxsu ) .* ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
b =  squeeze( sum( ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
c = a ./ b;
chlCarbonMxsuLat = squeeze( nanmean( chlCarbonMxsu ) );
chlCarbonMxsuLat( c < 0.5 ) = nan;

a = squeeze( sum( ~isnan( chlCarbonMxau ) .* ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
b =  squeeze( sum( ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
c = a ./ b;
chlCarbonMxauLat = squeeze( nanmean( chlCarbonMxau ) );
chlCarbonMxauLat( c < 0.5 ) = nan;

save( [PATH_OT, '/', FILENAME_OT], 'chlCarbonMxwi', 'chlCarbonMxsp', 'chlCarbonMxsu', 'chlCarbonMxau', 'chlCarbonMxwiLat', 'chlCarbonMxspLat', 'chlCarbonMxsuLat', 'chlCarbonMxauLat', 'mm', 'lat' )
disp( ['Save successfully: ', FILENAME_OT] )


%% 8. CHL
clear

PATH_IN = 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\biomes_chl_ecco';
PATH_OT = 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU';
FILENAME_OT = 'chlLat_season_ecco';

d = load( 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\biomes_npp_ecco\2003\npp_A EQU.mat' );
lon = d.lon;
lat = d.lat;
clear d;

d = load( 'm_topo.mat' );
topo = d.topo';

u1 = length( lon ); 
u2 = length( lat );
chl = nan( length( lon ), length( lat ), 20 * 12 );
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
            chl(I, J, (k - 1) * 12 + 1 : k * 12) = d.chl(i, :);
            basin_no(I, J) = basin;
        end
        basin = basin + 1;
    end
    basin = 0;
end

mm = repmat([1 : 12], 1, 20)';

ind1 = mm == 12 | mm == 1 | mm == 2; % Dec, Jan, Feb
chl1 = chl(:, :, ind1);
for i = 1 : 20
    chlwi(:, :, i) = nanmean( chl1(:, :, (i - 1) * 3 + 1 : (i - 1) * 3 + 3), 3);
end

ind2 = mm == 3 | mm == 4 | mm == 5;  % Mar, April, May
chl1 = chl(:, :, ind2);
for i = 1 : 20
    chlsp(:, :, i) = nanmean( chl1(:, :, (i - 1) * 3 + 1 : (i - 1) * 3 + 3), 3);
end

ind3 = mm == 6 | mm == 7 | mm == 8; % June, July, Aug
chl1 = chl(:, :, ind3);
for i = 1 : 20
    chlsu(:, :, i) = nanmean( chl1(:, :, (i - 1) * 3 + 1 : (i - 1) * 3 + 3), 3);
end

ind4 = mm == 9 | mm == 10 | mm == 11; % Sep, Oct, Nov
chl1 = chl(:, :, ind4);
for i = 1 : 20
    chlau(:, :, i) = nanmean( chl1(:, :, (i - 1) * 3 + 1 : (i - 1) * 3 + 3), 3);
end

idxN = lat > 0;
idxS = lat < 0;

chlwiS = chlwi(:, idxS, :);
chlwi(:, idxS, :) = chlsu(:, idxS, :);
chlsu(:, idxS, :) = chlwiS;

chlspS = chlsp(:, idxS, :);
chlsp(:, idxS, :) = chlau(:, idxS, :);
chlau(:, idxS, :) = chlspS;

a = squeeze( sum( ~isnan( chlwi ) .* ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
b =  squeeze( sum( ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
c = a ./ b;
chlwiLat = squeeze( nanmean( chlwi ) );
chlwiLat( c < 0.5 ) = nan;

a = squeeze( sum( ~isnan( chlsp ) .* ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
b =  squeeze( sum( ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
c = a ./ b;
chlspLat = squeeze( nanmean( chlsp ) );
chlspLat( c < 0.5 ) = nan;

a = squeeze( sum( ~isnan( chlsu ) .* ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
b =  squeeze( sum( ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
c = a ./ b;
chlsuLat = squeeze( nanmean( chlsu ) );
chlsuLat( c < 0.5 ) = nan;

a = squeeze( sum( ~isnan( chlau ) .* ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
b =  squeeze( sum( ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
c = a ./ b;
chlauLat = squeeze( nanmean( chlau ) );
chlauLat( c < 0.5 ) = nan;

save( [PATH_OT, '/', FILENAME_OT], 'chlwi', 'chlsp', 'chlsu', 'chlau', 'chlwiLat', 'chlspLat', 'chlsuLat', 'chlauLat', 'mm', 'lat' )
disp( ['Save successfully: ', FILENAME_OT] )

%% 9. PAR
clear

PATH_IN = 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\biomes_par_ecco';
PATH_OT = 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU';
FILENAME_OT = 'parLat_season_ecco';

d = load( 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\biomes_npp_ecco\2003\npp_A EQU.mat' );
lon = d.lon;
lat = d.lat;
clear d;

d = load( 'm_topo.mat' );
topo = d.topo';

u1 = length( lon ); 
u2 = length( lat );
par = nan( length( lon ), length( lat ), 20 * 12 );
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
            par(I, J, (k - 1) * 12 + 1 : k * 12) = d.par(i, :);
            basin_no(I, J) = basin;
        end
        basin = basin + 1;
    end
    basin = 0;
end

mm = repmat([1 : 12], 1, 20)';

ind1 = mm == 12 | mm == 1 | mm == 2; % Dec, Jan, Feb
par1 = par(:, :, ind1);
for i = 1 : 20
    parwi(:, :, i) = nanmean( par1(:, :, (i - 1) * 3 + 1 : (i - 1) * 3 + 3), 3);
end

ind2 = mm == 3 | mm == 4 | mm == 5;  % Mar, April, May
par1 = par(:, :, ind2);
for i = 1 : 20
    parsp(:, :, i) = nanmean( par1(:, :, (i - 1) * 3 + 1 : (i - 1) * 3 + 3), 3);
end

ind3 = mm == 6 | mm == 7 | mm == 8; % June, July, Aug
par1 = par(:, :, ind3);
for i = 1 : 20
    parsu(:, :, i) = nanmean( par1(:, :, (i - 1) * 3 + 1 : (i - 1) * 3 + 3), 3);
end

ind4 = mm == 9 | mm == 10 | mm == 11; % Sep, Oct, Nov
par1 = par(:, :, ind4);
for i = 1 : 20
    parau(:, :, i) = nanmean( par1(:, :, (i - 1) * 3 + 1 : (i - 1) * 3 + 3), 3);
end

idxN = lat > 0;
idxS = lat < 0;

parwiS = parwi(:, idxS, :);
parwi(:, idxS, :) = parsu(:, idxS, :);
parsu(:, idxS, :) = parwiS;

parspS = parsp(:, idxS, :);
parsp(:, idxS, :) = parau(:, idxS, :);
parau(:, idxS, :) = parspS;

a = squeeze( sum( ~isnan( parwi ) .* ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
b =  squeeze( sum( ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
c = a ./ b;
parwiLat = squeeze( nanmean( parwi ) );
parwiLat( c < 0.5 ) = nan;

a = squeeze( sum( ~isnan( parsp ) .* ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
b =  squeeze( sum( ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
c = a ./ b;
parspLat = squeeze( nanmean( parsp ) );
parspLat( c < 0.5 ) = nan;

a = squeeze( sum( ~isnan( parsu ) .* ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
b =  squeeze( sum( ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
c = a ./ b;
parsuLat = squeeze( nanmean( parsu ) );
parsuLat( c < 0.5 ) = nan;

a = squeeze( sum( ~isnan( parau ) .* ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
b =  squeeze( sum( ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
c = a ./ b;
parauLat = squeeze( nanmean( parau ) );
parauLat( c < 0.5 ) = nan;

save( [PATH_OT, '/', FILENAME_OT], 'parwi', 'parsp', 'parsu', 'parau', 'parwiLat', 'parspLat', 'parsuLat', 'parauLat', 'mm', 'lat' )
disp( ['Save successfully: ', FILENAME_OT] )