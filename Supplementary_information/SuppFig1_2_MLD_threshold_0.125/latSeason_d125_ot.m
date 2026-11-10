% Used for: seasonal latitudinal 0.125 MLD, mldLat_season_d125.mat (Supplementary Fig. 1)
% Original: Manuscript\MLD\manuscript\nature\NCC\R1\code for revision\mld_d125\latSeason_d125_ot.m (data paths updated to I:\Prj; computation unchanged)
% latSeason_d125_ot.m
%
% Stage 2 of the IDENTICAL two-stage pipeline used for mldLat_season
% (Code\mld_NEMO\latSeason_ot.m), applied to the 0.125 kg/m^3 MLD
% threshold. REPLACES an earlier version (2026-07-25) that pulled values
% directly from the full grid at reused 0.03 biome idx as a shortcut --
% per user request, this now reads from the freshly-regenerated
% biomes_1M_1dg\mld_d125\<year>\mld_<biome>.mat pre-split files (built by
% biomes_1M_d125_ot.m), exactly mirroring the 0.03 pipeline's two stages
% instead of skipping stage 1.
%
% Difference from Code\mld_NEMO\latSeason_ot.m (original untouched):
%   - PATH_IN: biomes_1M_1dg\mld_d125 (was biomes_1M_1dg\mld)
%   - lat/lon reference source: one_d125\mld\mld_2003.mat (was the
%     0.03-threshold npp_cafe biome file) -- same shared 1x1deg grid,
%     already verified elsewhere in this project to match exactly
%   - m_map path added explicitly for m_topo.mat
% Everything else (topo>-1000 depth mask, basin_no bookkeeping, DJF/MAM/
% JJA/SON seasonal averaging, hemisphere season-swap, 50% per-lat-band
% completeness threshold) is identical to the original.
%
% Output: mldLat_season_d125.mat (this folder).
%
% Haipeng Zhao
% 2026-07-25

clear; clc;
addpath('D:\Users\zhp19\OneDrive - Indiana State University\Research\MATLAB toolbox\m_map'); % m_topo.mat
addpath('D:\Users\zhp19\OneDrive - Indiana State University\Research\Projects\3D bio model\Code\mld_NEMO'); % m_topo.mat location fallback
OUT_DIR = fileparts(mfilename('fullpath'));

PATH_IN = 'I:\Prj\MLD_1M_NEMO\biomes_1M_1dg\mld_d125';

d = load('m_topo.mat');
topo = d.topo';

ref = load('I:\Prj\MLD_1M_NEMO\one_d125\mld\mld_2003.mat', 'lat','lon');
lon = ref.lon(:);
lat = ref.lat(:);

u1 = length( lon );
u2 = length( lat );
mld = nan( length( lon ), length( lat ), 20 * 12 );
basin_no = nan( length( lon ), length( lat ) );

datafile = dir( [PATH_IN, '\2003\', '*.mat'] );
filenum = length( datafile );
fprintf('Found %d biome files\n', filenum);

basin = 0;
for k = 1 : 20
    yyyy = num2str( k + 2002 );
    for j = 1 : filenum
        regfile = datafile(j);
        filename = regfile.name;
        disp( filename )

        d = load( fullfile( PATH_IN, yyyy, filename ) );
        for i = 1 : length( d.idx )
            [I, J] = ind2sub( [u1, u2], d.idx(i) );
            if ( topo(I, J) > -1000 )
                continue
            end
            mld(I, J, (k - 1) * 12 + 1 : k * 12) = d.mld(i, :);
            basin_no(I, J) = basin;
        end
        basin = basin + 1;
    end
    basin = 0;
end

mm = repmat( 1 : 12, 1, 20 )';

idxN = lat > 0;
idxS = lat < 0;

ind1 = mm == 12 | mm == 1 | mm == 2; % Dec, Jan, Feb
mld1 = mld(:, :, ind1);
for i = 1 : 20
    mldwi(:, :, i) = nanmean( mld1(:, :, (i - 1) * 3 + 1 : (i - 1) * 3 + 3), 3);
end

ind2 = mm == 3 | mm == 4 | mm == 5;  % Mar, April, May
mld1 = mld(:, :, ind2);
for i = 1 : 20
    mldsp(:, :, i) = nanmean( mld1(:, :, (i - 1) * 3 + 1 : (i - 1) * 3 + 3), 3);
end

ind3 = mm == 6 | mm == 7 | mm == 8; % June, July, Aug
mld1 = mld(:, :, ind3);
for i = 1 : 20
    mldsu(:, :, i) = nanmean( mld1(:, :, (i - 1) * 3 + 1 : (i - 1) * 3 + 3), 3);
end

ind4 = mm == 9 | mm == 10 | mm == 11; % Sep, Oct, Nov
mld1 = mld(:, :, ind4);
for i = 1 : 20
    mldau(:, :, i) = nanmean( mld1(:, :, (i - 1) * 3 + 1 : (i - 1) * 3 + 3), 3);
end

mldwiS = mldwi(:, idxS, :);
mldwi(:, idxS, :) = mldsu(:, idxS, :);
mldsu(:, idxS, :) = mldwiS;

mldspS = mldsp(:, idxS, :);
mldsp(:, idxS, :) = mldau(:, idxS, :);
mldau(:, idxS, :) = mldspS;

a = squeeze( sum( ~isnan( mldwi ) .* ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
b =  squeeze( sum( ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
c = a ./ b;
mldwiLat = squeeze( nanmean( mldwi ) );
mldwiLat( c < 0.5 ) = nan;

a = squeeze( sum( ~isnan( mldsp ) .* ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
b =  squeeze( sum( ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
c = a ./ b;
mldspLat = squeeze( nanmean( mldsp ) );
mldspLat( c < 0.5 ) = nan;

a = squeeze( sum( ~isnan( mldsu ) .* ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
b =  squeeze( sum( ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
c = a ./ b;
mldsuLat = squeeze( nanmean( mldsu ) );
mldsuLat( c < 0.5 ) = nan;

a = squeeze( sum( ~isnan( mldau ) .* ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
b =  squeeze( sum( ~isnan( repmat( basin_no, 1, 1, 20 ) ) ) );
c = a ./ b;
mldauLat = squeeze( nanmean( mldau ) );
mldauLat( c < 0.5 ) = nan;

save( fullfile(OUT_DIR, 'mldLat_season_d125'), 'mldwi', 'mldsp', 'mldsu', 'mldau', 'mldwiLat', 'mldspLat', 'mldsuLat', 'mldauLat', 'mm', 'lat' )
fprintf('Saved: %s\n', fullfile(OUT_DIR, 'mldLat_season_d125.mat'));
