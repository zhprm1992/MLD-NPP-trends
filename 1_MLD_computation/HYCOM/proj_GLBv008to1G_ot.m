% Used for: monthly MLD-HYCOM at 1 deg for 2003-2019, one_gat\mld\mld_YYYY.mat
% Original: Code\mld_HYCOM\proj_GLBv008to1G_ot.m (data paths updated to I:\Prj; computation unchanged)
% Project HYCOM GLBv008 to 1 degree
% Haipeng Zhao
% 09/12/2023
% 10/12/2023
% 9/26/2024
% 3/18/2025
% 6/21/2025

clear
clc

load( 'mldGLBv008to1GIdx.mat'  );
    
PATH = 'I:\Prj\MLD_1M_HYCOM_GLB\one12th'; % FOLDER NAME
PATH_OT = 'I:\Prj\MLD_HYCOM\one_gat\mld'; % FOLDER NAME

%% Target coordinates (ARGO JMSTEC)
fName = 'I:\Prj\MLD_ARGO\one\0.03sigmaTheta_ml_2003.nc\ml_20030101-20030131.nc';
lon = ncread( fName, 'LONGITUDE' );
lat = ncread( fName, 'LATITUDE' );
u1 = length( lon ); u2 = length( lat ); 

%%
[idx, ~, ic] = unique( gridIdx );
for i = 1 : 17
    yyyy = num2str( i + 2002 );
    mldG1deg = nan( u2, u1, 12 );
    
    file_mld = dir( [PATH, '\', yyyy, '\', '*.mat'] );
    for j = 1 : length( file_mld )
        token = split( file_mld(j).name, {'.', '_'} );
        mm = str2double( token{5} );

        load( [file_mld(j).folder, '\', file_mld(j).name], 'mld_1M' );
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
 


