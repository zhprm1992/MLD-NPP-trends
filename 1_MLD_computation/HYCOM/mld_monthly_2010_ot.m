% Used for: cluster (PACE) version of the yearly HYCOM MLD script; the year is set in the file
% Original: Code\mld_HYCOM\PACE\mld_monthly_2010_ot.m (data paths updated to I:\Prj; computation unchanged)
% Read HYCOM (reanalysis)
% Haipeng Zhao
% 6/6/2025

clear
clc

YYYY = 2010; % SET UP YEAR

DataFolderNm = ['./', num2str( YYYY )];

DataFolder = 'I:\Prj\MLD_1M_HYCOM\GLBv008';
DataFolder = ['./GLBv008_', num2str( YYYY )];

% ROOT_DAP_ULR = 'http://tds.hycom.org/thredds/dodsC';
ROOT_DAP_ULR = 'http://data.hycom.org/';

xmlURL = ['http://tds.hycom.org/thredds/catalog/datasets/GLBv0.08/expt_53.X/data/',...
    num2str( YYYY ), '/catalog.xml'];

% localFile = 'hycom_catalog.xml';
% websave( localFile, xmlURL );
doc = xmlread( xmlURL );

datasetNodes = doc.getElementsByTagName('dataset');
n = datasetNodes.getLength();

urlPaths = strings(0);
yyyy = [];
mm = [];
dd = [];
hh = [];
for i = 0 : n - 1
    node = datasetNodes.item(i);
    if node.hasAttribute( 'urlPath' )
        path = char( node.getAttribute( 'urlPath' ) );
        
        tokens = regexp( path, '\d{10}', 'match' );
        date_str = tokens{1};  % '2002123112'
        
        % Extract year, month, day
        yyyy(end + 1, 1)  = str2double( date_str(1 : 4));  % 
        mm(end + 1, 1) = str2double( date_str(5 : 6)); % 
        dd(end + 1, 1) = str2double( date_str(7 : 8)); % 
        hh(end + 1, 1) = str2double( path(end - 5 : end - 3) ); % 
        urlPaths(end + 1, 1) = path;
    end
end

num_lat = 3251;
num_lon = 4500;
num_dep = 40;

idlinkfull = "http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/";
ads1 = sprintf('%s%d', idlinkfull(1), yyyy );

ads5 = '?lat,lon';
ads_f5 = sprintf( '%s%s', ads1, ads5 );
lat = single( ncread( ads_f5, 'lat' ) );
lon = single( ncread( ads_f5, 'lon' ) );

ads3 = '?depth';
ads_f3 = sprintf( '%s%s', ads1, ads3 );
depth = single( ncread( ads_f3, 'depth' ) );

for i = 7 : 12
    n_days = eomday( YYYY, i );
    mld_dens = nan( num_lat, num_lon, n_days, 'single' );
    N200 = nan( num_lat, num_lon, n_days, 'single' );
    for j = 1 : n_days
        idx = find( (yyyy == YYYY) & (mm == i) & (dd == j) );
        path = urlPaths( idx );

        n_files = length( path );

        temp_1d = nan( num_lon, num_lat, num_dep, n_files, 'single' );
        sal_1d = nan( num_lon, num_lat, num_dep, n_files, 'single' );

        for k = 1 : n_files
            dataPATH = char( fullfile( ROOT_DAP_ULR, path(k) ) );
            filename = dataPATH( end - 36 : end );
            localFile = fullfile( DataFolder, filename );
            
            fprintf( 'Downloading file...' )
            cmd = sprintf( 'wget -q -O "%s" "%s"', localFile, dataPATH );
            status = system( cmd );

            if status == 0
                disp( dataPATH )
                temp_1d(:, :, :, k) = ncread( fullfile( DataFolder, filename ), 'water_temp' );
                sal_1d(:, :, :, k) = ncread( fullfile( DataFolder, filename ), 'salinity' );
            else
                warning('Fail to download file: %s\n', dataPATH);
            end    
        end
        delete( [DataFolder, '/*.nc'] );

        temp_1d = nanmean( temp_1d, 4 );
        sal_1d = nanmean( sal_1d, 4 );
 
        parfor m = 1 : num_lat
            for n = 1 : num_lon
                if ( isnan( temp_1d(n, m, 6) ) )
                    continue
                end
                temp_col = squeeze( temp_1d(n, m, :) );
                sal_col = squeeze( sal_1d(n, m, :) );

                PRES = gsw_p_from_z( -depth, lat(m) );
                SA = gsw_SA_from_SP( sal_col, PRES, lon(n), lat(m) );
                CT = gsw_CT_from_t( SA, temp_col, PRES );
                [N2, ~] = gsw_Nsquared( SA, CT, PRES, lat(m) );
              
                N200(m, n, j) = nanmean( N2((PRES >= 0) & (PRES <= 200)) );
                mld = findMLD_sigma003( CT, SA, PRES );
                if ~isempty( mld )
                    mld_dens(m, n, j) = mld;
                end
            end
        end
    end
    N200_1M = nanmean( N200, 3 );
    mld_1M = nanmean( mld_dens, 3 );
    filename_ot = ['mld_008_1M_', num2str( YYYY ), '_', sprintf( '%02d', i )];
    save( fullfile( DataFolderNm, filename_ot ), 'mld_1M', 'N200_1M', "mld_dens", 'N200', 'lat', 'lon', 'depth', '-v7.3' )
end





