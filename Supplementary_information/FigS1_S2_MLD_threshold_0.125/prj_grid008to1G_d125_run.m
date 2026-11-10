% Used for: regrid the 0.125 MLD to 1 deg (SI Figs. S1, S2)
% Original: Manuscript\MLD\manuscript\nature\NCC\R1\code for revision\mld_d125\prj_grid008to1G_d125_run.m (data paths updated to I:\Prj; computation unchanged)
% prj_grid008to1G_d125_run.m
%
% Runnable copy of prj_grid008to1G_d125_ot.m (original left untouched).
% Change: the target-coordinate source path was stale --
%   I:\Prj\MLD_ARGO\one\0.03sigmaTheta_ml_2003.nc\ml_20030101-20030131.nc
% no longer resolves (E:\ exists but has unrelated content; verified
% 2026-07-25). The same file now lives at:
%   I:\Prj\ml_20030101-20030131.nc
% Confirmed byte-identical LONGITUDE/LATITUDE to the file at the stale
% path's intended target AND to the shared 1x1deg grid used throughout
% this project (I:\Prj\MLD_1M_NEMO\one_d03\mld\mld_2003.mat, the
% corrected MLD source as of 2026-07-25 -- coordinates are identical
% between one\ and one_d03\, only the MLD values differ, so this
% coordinate check is unaffected by that bug) -- so this is a path
% fix only, not a coordinate change. mldG008to1G_idx.mat (same folder as
% this script) was double-checked too: gridIdx max=64800 = 360*180,
% consistent with this target grid.
%
% Haipeng Zhao
% 2026-07-25

clear
clc

load( 'mldG008to1G_idx.mat'  );

PATH = 'I:\Prj\MLD_1M_NEMO\one12th_d125'; % FOLDER NAME
PATH_OT = 'I:\Prj\MLD_1M_NEMO\one_d125\mld'; % FOLDER NAME
if ~exist(PATH_OT, 'dir'), mkdir(PATH_OT); end

%% Target coordinates (ARGO JMSTEC)
fName = 'I:\Prj\ml_20030101-20030131.nc';
lon = ncread( fName, 'LONGITUDE' );
lat = ncread( fName, 'LATITUDE' );
u1 = length( lon ); u2 = length( lat );

%%
[idx, ~, ic] = unique( gridIdx );
for i = 1 : 20   % SET UP MAX YEAR
    yyyy = num2str( i + 2002 );
    mldG1deg = nan( u2, u1, 12 );

    file_mld = dir( [PATH, '\', yyyy, '\', '*.nc'] );
    if numel(file_mld) < 12
        warning('%s: only %d/12 monthly files found, skipping year', yyyy, numel(file_mld));
        continue
    end
    for j = 1 : 12
        mld_1M = ncread( [file_mld(j).folder, '\', file_mld(j).name], 'mld_d125' );
        mld_1M = mld_1M';
        mld = nan( u2, u1 );
        mld(idx) = splitapply( @(x) mean( x, 'omitnan' ), mld_1M(:), ic );
        mldG1deg(:, :, j) = mld;
    end

    %%
    % Output variables
    filename_ot = ['mld_', yyyy]; % FILE NAME
    save(  [PATH_OT, '\', filename_ot],  'mldG1deg', 'lon', 'lat', '-v7.3'  ) % Grid = [90, 180]
    fprintf('Saved: %s\n', [PATH_OT, '\', filename_ot]);
end
