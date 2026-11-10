% Used for: monthly MLD-NEMO from 1/12 deg to 1 deg (mean of the 1/12 deg cells), one_d03\mld\mld_YYYY.mat; set for MLD (variable mld_d03), the same script was also run for N2 (bvf_200)
% Original: Code\mld_NEMO\prj_grid008to1G_ot_v2.m (data paths updated to I:\Prj; computation unchanged)
% Project 1/6  degree to 1 degree. This version uses MLD density 0.03 km/m3
% Haipeng Zhao
% 09/12/2023
% 10/12/2023
% 9/26/2024
% 3/18/2025

clear
clc

load( 'mldG008to1G_idx.mat'  );
    
PATH = 'I:\Prj\mld_1M_NEMO\one12th_d03'; % FOLDER NAME
PATH_OT = 'I:\Prj\mld_1M_NEMO\one_d03\mld'; % FOLDER NAME

%% Target coordinates (ARGO JMSTEC)
fName = 'I:\Prj\MLD_ARGO\one\0.03sigmaTheta_ml_2003.nc\ml_20030101-20030131.nc';
lon = ncread( fName, 'LONGITUDE' );
lat = ncread( fName, 'LATITUDE' );
u1 = length( lon ); u2 = length( lat ); 

%%
[idx, ~, ic] = unique( gridIdx );
for i = 1 : 20   % SET UP MAX YEAR
    yyyy = num2str( i + 2002 );
    mldG1deg = nan( u2, u1, 12 );
    
    file_mld = dir( [PATH, '\', yyyy, '\', '*.nc'] );
    for j = 1 : 12
        mld_1M = ncread( [file_mld(j).folder, '\', file_mld(j).name], 'mld_d03' );
        mld_1M = mld_1M';
        mld = nan( u2, u1 );   
        mld(idx) = splitapply( @(x) mean( x, 'omitnan' ), mld_1M(:), ic );
        mldG1deg(:, :, j) = mld;
    end

    %%
    % Output variables
    filename_ot = ['mld_', yyyy]; % FILE NAME
    save(  [PATH_OT, '\', filename_ot],  'mldG1deg', 'lon', 'lat', '-v7.3'  ) % Grid = [90, 180]
end
 
 


