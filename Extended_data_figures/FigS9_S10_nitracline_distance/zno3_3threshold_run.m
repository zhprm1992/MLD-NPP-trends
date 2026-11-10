% Used for: WOA23 nitracline depth (3 umol/kg), zno3_woa23_3muM_nofix.mat (Figs. S9, S10)
% Original: Code\npp_cbpm\zno3_3threshold_run.m (data paths updated to I:\Prj; computation unchanged)
% zno3_3threshold_run.m
%
% Runnable copy of zno3_3threshold_ot.m (original untouched). Removes the
% one line that substitutes a fixed 333m value south of 30S (the Arteaga
% et al. 2022 ferricline convention) -- this version computes the actual
% NO3=3 umol/kg threshold depth everywhere, including south of 30S, so
% the full field can be mapped without an artificial flat patch. This is
% NOT the field used by the CbPM NPP generator (which needs the Arteaga
% substitution for the ferricline reasoning south of 30S) -- it exists
% purely to visualize what the actual nitracline looks like everywhere.
%
% Output: zno3_woa23_3muM_nofix.mat (I:\Prj\WOA23\nitracline\, same folder
% as the original) -- same fields (zno3, lat, lon), just without the
% south-of-30S substitution.
%
% Haipeng Zhao
% 2026-07-28

clear
clc

PATH = 'I:\Prj\WOA23\nitrate';
PATH_OT = 'I:\Prj\WOA23\nitracline';

datafile = dir( fullfile( PATH, '*nc') );

for i = 1 : 12
    zno_depth = nan( 180 , 360, 'single' ); % depth of the nitracline

    zno = ncread( fullfile( datafile(i).folder, datafile(i).name ), 'n_an' ); % mu mol/kg
    zno = permute( zno, [2, 1, 3] );
    depth = ncread( fullfile( datafile(i).folder, datafile(i).name ), 'depth' );
    lat = ncread( fullfile( datafile(i).folder, datafile(i).name ), 'lat' );
    lon = ncread( fullfile( datafile(i).folder, datafile(i).name ), 'lon' );

    for m = 1 : 180
        for n = 1 : 360
            no3_prof = squeeze( zno(m, n, :) );
            if all( isnan( no3_prof ) )
                continue
            end
            zno_depth(m, n) = get_nitracline( no3_prof, depth, 3 );
        end
    end
    zno3(:, :, i) = zno_depth;
    fprintf('month %d done\n', i);
end
% No south-of-30S substitution here (unlike zno3_3threshold_ot.m) --
% this is the actual computed nitracline depth everywhere.

filename_ot = 'zno3_woa23_3muM_nofix';
save( fullfile( PATH_OT, filename_ot ), 'zno3', 'lat', 'lon', '-v7.3' )
disp('Saved successfully: zno3 (no south-of-30S fix)');
