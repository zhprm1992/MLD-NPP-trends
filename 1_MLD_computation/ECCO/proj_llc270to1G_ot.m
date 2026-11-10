% Used for: monthly MLD-ECCO at 1 deg, mld_one\mld_YYYY.mat (copied to I:\Prj\MLD_ECCO\one)
% Original: Code\mld_ECCO\Compute_MLD\proj_llc270to1G_ot.m (data paths updated to I:\Prj; computation unchanged)
% Project HYCOM GLBv008 to 1 degree
% Haipeng Zhao
% 09/12/2023
% 10/12/2023
% 9/26/2024
% 3/18/2025
% 6/21/2025

clear
clc

load( 'mldLLC270to2G_idx.mat'  );
    
PATH = 'I:\ECCO-Darwin_extension\daily\mld_003\Monthly'; % FOLDER NAME
PATH_OT = 'I:\ECCO-Darwin_extension\daily\mld_003\mld_one'; % FOLDER NAME

numFaces = 13;
nx = 270;
ny = 270 .* numFaces;
nz = 50;

%% Target coordinates (ARGO JMSTEC)
fName = 'I:\Prj\MLD_ARGO\one\0.03sigmaTheta_ml_2003.nc\ml_20030101-20030131.nc';
lon = ncread( fName, 'LONGITUDE' );
lat = ncread( fName, 'LATITUDE' );
u1 = length( lon ); u2 = length( lat ); 

%%
[idx, ~, ic] = unique( gridIdx );
for i = 1 : 20
    yyyy = num2str( i + 2002 );
    mldG1deg = nan( u2, u1, 12 );
    
    fileList = dir( [PATH, '\', yyyy, '\', '*.data'] );
    for j = 1 : length( fileList )
        token = split( fileList(j).name, {'.', '_'} );
        mm = str2double( token{3} );

        mld_1M = readbin( fullfile( fileList(j).folder, fileList(j).name ), [nx ny], 1, 'real*4' );
        mld_1M(mld_1M < 0) = nan;
        mld = nan( u2, u1 );   
        mld(idx) = splitapply( @(x) mean( x, 'omitnan' ), mld_1M(:), ic );
        mldG1deg(:, :, mm) = mld;
    end

    %%
    % Output variables
    filename_ot = ['mld_', yyyy]; % FILE NAME
    save(  [PATH_OT, '\', filename_ot],  'mldG1deg', 'lon', 'lat', '-v7.3'  ) % Grid = [90, 180]
end
 


