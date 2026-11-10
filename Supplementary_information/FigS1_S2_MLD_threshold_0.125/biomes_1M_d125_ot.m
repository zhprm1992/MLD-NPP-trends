% Used for: 0.125 MLD by biome
% Original: Manuscript\MLD\manuscript\nature\NCC\R1\code for revision\mld_d125\biomes_1M_d125_ot.m (data paths updated to I:\Prj; computation unchanged)
% biomes_1M_d125_ot.m
%
% Stage 1 of the IDENTICAL two-stage pipeline used to build mldLat_season
% (Code\mld_NEMO\biomes_1M_ot.m + latSeason_ot.m), applied to the 0.125
% kg/m^3 MLD threshold instead of the manuscript's 0.03 threshold. Splits
% the full 1x1deg grid into per-biome-per-year files by re-running the
% SAME inpolygon biome-split as the original (fresh idx, not reused from
% the 0.03 pre-split files), so the d125 latitudinal-average pipeline is
% byte-for-byte the same procedure as d03, not a shortcut. Per user
% request 2026-07-25 -- an earlier version (latSeason_d125_ot.m, since
% replaced) pulled values directly from the full grid at reused 0.03 idx
% instead of regenerating this intermediate stage; mathematically expected
% to be equivalent but not verified byte-for-byte, and this project has
% already been burned once by an unverified "should be equivalent"
% shortcut (the one\ vs one_d03\ MLD source bug) -- so this stage is now
% reproduced explicitly instead of assumed.
%
% Difference from Code\mld_NEMO\biomes_1M_ot.m (original untouched):
%   - PATH_IN:  one_d125\mld  (was one_d03\mld)
%   - PATH_OT:  biomes_1M_1dg\mld_d125  (was biomes_1M_1dg\mld -- new
%     sibling folder, does NOT overwrite the 0.03 pre-split files)
%   - lat/lon for the inpolygon grid read directly from one_d125's own
%     mld_2003.mat (already verified elsewhere in this project to match
%     the shared 1x1deg grid exactly) instead of the stale
%     I:\Prj\MLD_ARGO\one\0.03sigmaTheta_ml_2003.nc\... path
%   - SEABIOMES shapefile repointed to its current location (D:\...\
%     DINEOF application\biomes\meanBiom_FAY\meanBiom_FAY.shp -- E:\ no
%     longer has this file, same drive-letter drift seen elsewhere)
%   - Filenames constructed directly (mld_YYYY.mat) instead of relying on
%     dir() listing order
% Everything else (inpolygon biome-split logic, idx orientation, mld
% extraction) is identical to the original.
%
% NOTE (2026-07-25): the original biomes_1M_ot.m on disk builds
% fileName = biomename (no prefix), which would save as e.g. "A EQU.mat".
% But the actual 0.03 pre-split files on disk are named "mld_A EQU.mat",
% and every downstream script in this project (regress_biome_mld_npp_
% trends.m, latSeason_ot.m via dir()'s *.mat glob, etc.) is agnostic to
% or expects that mld_ prefix. The script-as-currently-saved is therefore
% NOT byte-identical to whatever actually produced the existing files
% (it may have been edited after the fact) -- this does not affect the
% numeric content (inpolygon/idx/mld extraction is unambiguous), only the
% output filename. This script uses fileName = ['mld_', biomename] to
% match the real, working convention every other script in this project
% relies on.
%
% Haipeng Zhao
% 2026-07-25

clear
clc

PATH_IN = 'I:\Prj\MLD_1M_NEMO\one_d125\mld';
PATH_OT = 'I:\Prj\MLD_1M_NEMO\biomes_1M_1dg\mld_d125';
SEABIOMES = shaperead( ['D:\Users\zhp19\OneDrive - Indiana State University\Research\Projects'...
    '\', 'DINEOF application\biomes\meanBiom_FAY',...
    '\', 'meanBiom_FAY.shp'] );

addpath('D:\Users\zhp19\OneDrive - Indiana State University\Research\Projects\3D bio model\Code\mld_NEMO'); % nameConvert.m

ref = load( fullfile( PATH_IN, 'mld_2003.mat' ), 'lat', 'lon' );
lat = ref.lat;
lon = ref.lon;

for i = 1 : 20
    yyyy = num2str( i + 2002 );

    load( fullfile( PATH_IN, ['mld_', yyyy, '.mat'] ) )
    [lat_g, lon_g] = meshgrid( lat, lon );
    lon_g = single( lon_g(:) );
    lat_g = single( lat_g(:) );

    for k = 1 : length( SEABIOMES )
        seaBd = SEABIOMES(k);
        idx = inpolygon( lon_g, lat_g, seaBd.X, seaBd.Y ); % find pixels inside of the region
        idx = find( idx );
        pixelnum = length( idx );
        area = seaBd.Label;
        biomename = nameConvert( area );
        mld = zeros( pixelnum, 12, 'single' );
        for j = 1 : 12
            mld1 = mldG1deg(:, :, j);
            mld1 = mld1';
            mld(:, j) = mld1(idx);
        end

        fileName = ['mld_', biomename];
        path_ot = fullfile( PATH_OT, yyyy );
        if ~exist( path_ot, 'dir' )
            mkdir( path_ot );
        end
        save( [path_ot, '\', fileName], 'mld', 'lat', 'lon', 'idx', '-v7.3' );
    end
    fprintf('%s done\n', yyyy);
end
