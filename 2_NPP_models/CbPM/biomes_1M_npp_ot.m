% Used for: CbPM output split into biomes (NPP, growth rate, carbon, ...), biomes_npp_<product>; set here for ECCO, run once per MLD product
% Original: Code\npp_cbpm\biomes_1M_npp_ot.m (data paths updated to I:\Prj; computation unchanged)
% SechlCarbonate global NPlgP ilgo biomes
% Haipeng Zhao
% 9/9/2023
% 9/26/2023
% 10/11/2023 Remove values more than  three standard deviation from mean
% 10/17/2023 Extended to process SST dataset
% 09/25/2024 Applied to POC dataset
% 10/4/2024 Applied to NPP dataset

clc

%% 1. NPP
clear

PATH_IN = 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\grid_npp_ecco'; % CHECK INPUT FOLDER
PATH_OT = 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\biomes_npp_ecco'; % CHECK OUTPUT FOLDER
SEABIOMES = shaperead( ['..\..\..'...
    '\', 'DINEOF application\biomes\meanBiom_FAY',...
    '\', 'meanBiom_FAY.shp'] );

fName = 'I:\Prj\ml_20030101-20030131.nc';
lon = ncread( fName, 'LONGITUDE' );
lat = ncread( fName, 'LATITUDE' );

datafile = dir( fullfile( PATH_IN, '*.mat') );
for i = 1 : 20
    yyyy = num2str( i + 2002 );

    load( fullfile( datafile(i).folder, datafile(i).name ) )
    [lat_g, lon_g] = meshgrid( lat, lon );
    lon_g = single( lon_g(:) );
    lat_g = single( lat_g(:) );

    if ~exist( 'nppG1dg', 'var' )
        nppG1dg = npp;   % CHECK VARIABLE NAME
    end

    for k = 1 : length( SEABIOMES )
        seaBd = SEABIOMES(k);
        idx = inpolygon( lon_g, lat_g, seaBd.X, seaBd.Y ); % find pixels inside of the region   
        idx = find( idx );
        pixelnum = length( idx );
        area = seaBd.Label;
        biomename = nameConvert( area );
        npp = nan( pixelnum, 12, 'single' );
        for j = 1 : 12
            var1 = nppG1dg(:, :, j);
            var1 = var1';
            tmp = var1(idx);
            [B, TF] = rmoutliers( tmp );
            tmp(TF) = nan;
            npp(:, j) = tmp;            
        end
        fileName = ['npp_', biomename];  % CHECK FILENAME
        path_ot = fullfile( PATH_OT, yyyy );
        if ~exist( path_ot, 'dir' )
            mkdir( path_ot );
        end
        save( [path_ot, '\', fileName], 'npp', 'lat', 'lon', 'idx', '-v7.3' );
    end
end
disp('Saved successfully: NPP');


%% 2. MU
clear

PATH_IN = 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\grid_npp_ecco'; % CHECK INPUT FOLDER
PATH_OT = 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\biomes_mu_ecco'; % CHECK OUTPUT FOLDER
SEABIOMES = shaperead( ['..\..\..'...
    '\', 'DINEOF application\biomes\meanBiom_FAY',...
    '\', 'meanBiom_FAY.shp'] );

fName = 'I:\Prj\ml_20030101-20030131.nc';
lon = ncread( fName, 'LONGITUDE' );
lat = ncread( fName, 'LATITUDE' );

datafile = dir( fullfile( PATH_IN, '*.mat') );
for i = 1 : 20
    yyyy = num2str( i + 2002 );

    load( fullfile( datafile(i).folder, datafile(i).name ) )
    [lat_g, lon_g] = meshgrid( lat, lon );
    lon_g = single( lon_g(:) );
    lat_g = single( lat_g(:) );

    if ~exist( 'muG1dg', 'var' )
        muG1dg = mu;   % CHECK VARIABLE NAME
    end

    for k = 1 : length( SEABIOMES )
        seaBd = SEABIOMES(k);
        idx = inpolygon( lon_g, lat_g, seaBd.X, seaBd.Y ); % find pixels inside of the region   
        idx = find( idx );
        pixelnum = length( idx );
        area = seaBd.Label;
        biomename = nameConvert( area );
        mu = zeros( pixelnum, 12, 'single' );
        for j = 1 : 12
            var1 = muG1dg(:, :, j);
            var1 = var1';
            tmp = var1(idx);
            [B, TF] = rmoutliers( tmp );
            tmp(TF) = nan;
            mu(:, j) = tmp;               
        end

        fileName = ['mu_', biomename];  % CHECK FILENAME
        path_ot = fullfile( PATH_OT, yyyy );
        if ~exist( path_ot, 'dir' )
            mkdir( path_ot );
        end
        save( [path_ot, '\', fileName], 'mu', 'lat', 'lon', 'idx', '-v7.3' );
    end
end
disp('Saved successfully: MU');

%% 3. CARBON
clear

PATH_IN = 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\grid_npp_ecco'; % CHECK INPUT FOLDER
PATH_OT = 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\biomes_carbon_ecco'; % CHECK OUTPUT FOLDER
SEABIOMES = shaperead( ['..\..\..'...
    '\', 'DINEOF application\biomes\meanBiom_FAY',...
    '\', 'meanBiom_FAY.shp'] );

fName = 'I:\Prj\ml_20030101-20030131.nc';
lon = ncread( fName, 'LONGITUDE' );
lat = ncread( fName, 'LATITUDE' );

datafile = dir( fullfile( PATH_IN, '*.mat') );
for i = 1 : 20
    yyyy = num2str( i + 2002 );

    load( fullfile( datafile(i).folder, datafile(i).name ) )
    [lat_g, lon_g] = meshgrid( lat, lon );
    lon_g = single( lon_g(:) );
    lat_g = single( lat_g(:) );

    if ~exist( 'carbonG1dg', 'var' )
        carbonG1dg = carbon;   % CHECK VARIABLE NAME
    end

    for k = 1 : length( SEABIOMES )
        seaBd = SEABIOMES(k);
        idx = inpolygon( lon_g, lat_g, seaBd.X, seaBd.Y ); % find pixels inside of the region   
        idx = find( idx );
        pixelnum = length( idx );
        area = seaBd.Label;
        biomename = nameConvert( area );
        carbon = zeros( pixelnum, 12, 'single' );
        for j = 1 : 12
            var1 = carbonG1dg(:, :, j);
            var1 = var1';
            tmp = var1(idx);
            [B, TF] = rmoutliers( tmp );
            tmp(TF) = nan;
            carbon(:, j) = tmp;            
        end

        fileName = ['carbon_', biomename];  % CHECK FILENAME
        path_ot = fullfile( PATH_OT, yyyy );
        if ~exist( path_ot, 'dir' )
            mkdir( path_ot );
        end
        save( [path_ot, '\', fileName], 'carbon', 'lat', 'lon', 'idx', '-v7.3' );
        
    end
end
disp('Saved successfully: Carbon');

%% 4. IgFUN
clear

PATH_IN = 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\grid_npp_ecco'; % CHECK INPUT FOLDER
PATH_OT = 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\biomes_IgFun_ecco'; % CHECK OUTPUT FOLDER
SEABIOMES = shaperead( ['..\..\..'...
    '\', 'DINEOF application\biomes\meanBiom_FAY',...
    '\', 'meanBiom_FAY.shp'] );

fName = 'I:\Prj\ml_20030101-20030131.nc';
lon = ncread( fName, 'LONGITUDE' );
lat = ncread( fName, 'LATITUDE' );

datafile = dir( fullfile( PATH_IN, '*.mat') );
for i = 1 : 20
    yyyy = num2str( i + 2002 );

    load( fullfile( datafile(i).folder, datafile(i).name ) )
    [lat_g, lon_g] = meshgrid( lat, lon );
    lon_g = single( lon_g(:) );
    lat_g = single( lat_g(:) );

    if ~exist( 'IgFUNG1dg', 'var' )
        IgFUNG1dg = IgFun;   % CHECK VARIABLE NAME
    end

    for k = 1 : length( SEABIOMES )
        seaBd = SEABIOMES(k);
        idx = inpolygon( lon_g, lat_g, seaBd.X, seaBd.Y ); % find pixels inside of the region   
        idx = find( idx );
        pixelnum = length( idx );
        area = seaBd.Label;
        biomename = nameConvert( area );
        IgFun = zeros( pixelnum, 12, 'single' );
        for j = 1 : 12
            var1 = IgFUNG1dg(:, :, j);
            var1 = var1';
            tmp = var1(idx);
            [B, TF] = rmoutliers( tmp );
            tmp(TF) = nan;
            IgFun(:, j) = tmp;            
        end

        fileName = ['IgFun_', biomename];  % CHECK FILENAME
        path_ot = fullfile( PATH_OT, yyyy );
        if ~exist( path_ot, 'dir' )
            mkdir( path_ot );
        end
        save( [path_ot, '\', fileName], 'IgFun', 'lat', 'lon', 'idx', '-v7.3' );
    end
end
disp('Saved successfully: IgFun');

%% 5. ntFUN
clear

PATH_IN = 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\grid_npp_ecco'; % CHECK INPUT FOLDER
PATH_OT = 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\biomes_ntFun_ecco'; % CHECK OUTPUT FOLDER
SEABIOMES = shaperead( ['..\..\..'...
    '\', 'DINEOF application\biomes\meanBiom_FAY',...
    '\', 'meanBiom_FAY.shp'] );

fName = 'I:\Prj\ml_20030101-20030131.nc';
lon = ncread( fName, 'LONGITUDE' );
lat = ncread( fName, 'LATITUDE' );

datafile = dir( fullfile( PATH_IN, '*.mat') );
for i = 1 : 20
    yyyy = num2str( i + 2002 );

    load( fullfile( datafile(i).folder, datafile(i).name ) )
    [lat_g, lon_g] = meshgrid( lat, lon );
    lon_g = single( lon_g(:) );
    lat_g = single( lat_g(:) );

    if ~exist( 'ntFUNG1dg', 'var' )
        ntFUNG1dg = ntFun;   % CHECK VARIABLE NAME
    end

    for k = 1 : length( SEABIOMES )
        seaBd = SEABIOMES(k);
        idx = inpolygon( lon_g, lat_g, seaBd.X, seaBd.Y ); % find pixels inside of the region   
        idx = find( idx );
        pixelnum = length( idx );
        area = seaBd.Label;
        biomename = nameConvert( area );
        ntFun = zeros( pixelnum, 12, 'single' );
        for j = 1 : 12
            var1 = ntFUNG1dg(:, :, j);
            var1 = var1';
            tmp = var1(idx);
            [B, TF] = rmoutliers( tmp );
            tmp(TF) = nan;
            ntFun(:, j) = tmp;            
        end

        fileName = ['ntFun_', biomename];  % CHECK FILENAME
        path_ot = fullfile( PATH_OT, yyyy );
        if ~exist( path_ot, 'dir' )
            mkdir( path_ot );
        end
        save( [path_ot, '\', fileName], 'ntFun', 'lat', 'lon', 'idx', '-v7.3' );
    end
end
disp('Saved successfully: ntFun');

%% 6. chlCarbon
clear

PATH_IN = 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\grid_npp_ecco'; % CHECK INPUT FOLDER
PATH_OT = 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\biomes_chlCarbon_ecco'; % CHECK OUTPUT FOLDER
SEABIOMES = shaperead( ['..\..\..'...
    '\', 'DINEOF application\biomes\meanBiom_FAY',...
    '\', 'meanBiom_FAY.shp'] );

fName = 'I:\Prj\ml_20030101-20030131.nc';
lon = ncread( fName, 'LONGITUDE' );
lat = ncread( fName, 'LATITUDE' );

datafile = dir( fullfile( PATH_IN, '*.mat') );
for i = 1 : 20
    yyyy = num2str( i + 2002 );

    load( fullfile( datafile(i).folder, datafile(i).name ) )
    [lat_g, lon_g] = meshgrid( lat, lon );
    lon_g = single( lon_g(:) );
    lat_g = single( lat_g(:) );

    if ~exist( 'chlCarbonG1dg', 'var' )
        chlCarbonG1dg = chlCarbon;   % CHECK VARIABLE NAME
    end

    for k = 1 : length( SEABIOMES )
        seaBd = SEABIOMES(k);
        idx = inpolygon( lon_g, lat_g, seaBd.X, seaBd.Y ); % find pixels inside of the region   
        idx = find( idx );
        pixelnum = length( idx );
        area = seaBd.Label;
        biomename = nameConvert( area );
        chlCarbon = zeros( pixelnum, 12, 'single' );
        for j = 1 : 12
            var1 = chlCarbonG1dg(:, :, j);
            var1 = var1';
            tmp = var1(idx);
            [B, TF] = rmoutliers( tmp );
            tmp(TF) = nan;
            chlCarbon(:, j) = tmp;            
        end

        fileName = ['chlCarbon_', biomename];  % CHECK FILENAME
        path_ot = fullfile( PATH_OT, yyyy );
        if ~exist( path_ot, 'dir' )
            mkdir( path_ot );
        end
        save( [path_ot, '\', fileName], 'chlCarbon', 'lat', 'lon', 'idx', '-v7.3' );
    end
end
disp('Saved successfully: chlCarbon');

%% 7. chlCarbonMax
clear

PATH_IN = 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\grid_npp_ecco'; % CHECK INPUT FOLDER
PATH_OT = 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\biomes_chlCarbonMx_ecco'; % CHECK OUTPUT FOLDER
SEABIOMES = shaperead( ['..\..\..'...
    '\', 'DINEOF application\biomes\meanBiom_FAY',...
    '\', 'meanBiom_FAY.shp'] );

fName = 'I:\Prj\ml_20030101-20030131.nc';
lon = ncread( fName, 'LONGITUDE' );
lat = ncread( fName, 'LATITUDE' );

datafile = dir( fullfile( PATH_IN, '*.mat') );
for i = 1 : 20
    yyyy = num2str( i + 2002 );

    load( fullfile( datafile(i).folder, datafile(i).name ) )
    [lat_g, lon_g] = meshgrid( lat, lon );
    lon_g = single( lon_g(:) );
    lat_g = single( lat_g(:) );

    if ~exist( 'chlCarbonMxG1dg', 'var' )
        chlCarbonMxG1dg = chlCarbonMx;   % CHECK VARIABLE NAME
    end

    for k = 1 : length( SEABIOMES )
        seaBd = SEABIOMES(k);
        idx = inpolygon( lon_g, lat_g, seaBd.X, seaBd.Y ); % find pixels inside of the region   
        idx = find( idx );
        pixelnum = length( idx );
        area = seaBd.Label;
        biomename = nameConvert( area );
        chlCarbonMx = zeros( pixelnum, 12, 'single' );
        for j = 1 : 12
            var1 = chlCarbonMxG1dg(:, :, j);
            var1 = var1';
            tmp = var1(idx);
            [B, TF] = rmoutliers( tmp );
            tmp(TF) = nan;
            chlCarbonMx(:, j) = tmp;               
        end

        fileName = ['chlCarbonMx_', biomename];  % CHECK FILENAME
        path_ot = fullfile( PATH_OT, yyyy );
        if ~exist( path_ot, 'dir' )
            mkdir( path_ot );
        end
        save( [path_ot, '\', fileName], 'chlCarbonMx', 'lat', 'lon', 'idx', '-v7.3' );
    end
end
disp('Saved successfully: chlCarbonMx');


%% 8. chl
clear

PATH_IN = 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\grid_npp_ecco'; % CHECK INPUT FOLDER
PATH_OT = 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\biomes_chl_ecco'; % CHECK OUTPUT FOLDER
SEABIOMES = shaperead( ['..\..\..'...
    '\', 'DINEOF application\biomes\meanBiom_FAY',...
    '\', 'meanBiom_FAY.shp'] );

fName = 'I:\Prj\ml_20030101-20030131.nc';
lon = ncread( fName, 'LONGITUDE' );
lat = ncread( fName, 'LATITUDE' );

datafile = dir( fullfile( PATH_IN, '*.mat') );
for i = 1 : 20
    yyyy = num2str( i + 2002 );

    load( fullfile( datafile(i).folder, datafile(i).name ) )
    [lat_g, lon_g] = meshgrid( lat, lon );
    lon_g = single( lon_g(:) );
    lat_g = single( lat_g(:) );

    if ~exist( 'chlG1dg', 'var' )
        chlG1dg = chl;   % CHECK VARIABLE NAME
    end

    for k = 1 : length( SEABIOMES )
        seaBd = SEABIOMES(k);
        idx = inpolygon( lon_g, lat_g, seaBd.X, seaBd.Y ); % find pixels inside of the region   
        idx = find( idx );
        pixelnum = length( idx );
        area = seaBd.Label;
        biomename = nameConvert( area );
        chl = zeros( pixelnum, 12, 'single' );
        for j = 1 : 12
            var1 = chlG1dg(:, :, j);
            var1 = var1';
            tmp = var1(idx);
            [B, TF] = rmoutliers( tmp );
            tmp(TF) = nan;
            chl(:, j) = tmp;               
        end

        fileName = ['chl_', biomename];  % CHECK FILENAME
        path_ot = fullfile( PATH_OT, yyyy );
        if ~exist( path_ot, 'dir' )
            mkdir( path_ot );
        end
        save( [path_ot, '\', fileName], 'chl', 'lat', 'lon', 'idx', '-v7.3' );
    end
end
disp('Saved successfully: chl');


%% 9. par
clear

PATH_IN = 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\grid_npp_ecco'; % CHECK INPUT FOLDER
PATH_OT = 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\biomes_par_ecco'; % CHECK OUTPUT FOLDER
SEABIOMES = shaperead( ['..\..\..'...
    '\', 'DINEOF application\biomes\meanBiom_FAY',...
    '\', 'meanBiom_FAY.shp'] );

fName = 'I:\Prj\ml_20030101-20030131.nc';
lon = ncread( fName, 'LONGITUDE' );
lat = ncread( fName, 'LATITUDE' );

datafile = dir( fullfile( PATH_IN, '*.mat') );
for i = 1 : 20
    yyyy = num2str( i + 2002 );

    load( fullfile( datafile(i).folder, datafile(i).name ) )
    [lat_g, lon_g] = meshgrid( lat, lon );
    lon_g = single( lon_g(:) );
    lat_g = single( lat_g(:) );

    if ~exist( 'parG1dg', 'var' )
        parG1dg = par;   % CHECK VARIABLE NAME
    end

    for k = 1 : length( SEABIOMES )
        seaBd = SEABIOMES(k);
        idx = inpolygon( lon_g, lat_g, seaBd.X, seaBd.Y ); % find pixels inside of the region   
        idx = find( idx );
        pixelnum = length( idx );
        area = seaBd.Label;
        biomename = nameConvert( area );
        par = zeros( pixelnum, 12, 'single' );
        for j = 1 : 12
            var1 = parG1dg(:, :, j);
            var1 = var1';
            tmp = var1(idx);
            [B, TF] = rmoutliers( tmp );
            tmp(TF) = nan;
            par(:, j) = tmp;               
        end

        fileName = ['par_', biomename];  % CHECK FILENAME
        path_ot = fullfile( PATH_OT, yyyy );
        if ~exist( path_ot, 'dir' )
            mkdir( path_ot );
        end
        save( [path_ot, '\', fileName], 'par', 'lat', 'lon', 'idx', '-v7.3' );
    end
end
disp('Saved successfully: par');
