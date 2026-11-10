% Used for: monthly latitudinal 0.125 MLD, mldLat_monthly_d125.mat (SI Fig. S2)
% Original: Manuscript\MLD\manuscript\nature\NCC\R1\code for revision\mld_d125\monthlyLat_d125_ot.m (data paths updated to I:\Prj; computation unchanged)
% monthlyLat_d125_ot.m
%
% 0.125 kg/m^3 threshold equivalent of Code\mld_NEMO\monthlyLat_ot.m (the
% script that built mldLat_monthly_nemo.mat, the actual input consumed by
% the manuscript's real Figure 4 pipeline: fs6_heatmap_monthlyNPPtrends_ot.m
% -> NPPtrends.mat -> senMLD2_ot.m -> f4_plot_sensitivity_NPP_MLD.m).
% Reads from the freshly-regenerated biomes_1M_1dg\mld_d125\<year>\
% mld_<biome>.mat pre-split files (built by biomes_1M_d125_ot.m -- the
% same fresh inpolygon biome-split used for mldLat_season_d125.mat), same
% topo>-1000 depth mask, same basin bookkeeping structure as the original
% -- differs only in that it keeps all 240 individual year-months
% zonally-averaged (no 3-month seasonal averaging, no hemisphere swap),
% matching monthlyLat_ot.m exactly.
%
% Output: mldLat_monthly_d125.mat (this folder) -- mldLat (nlon-collapsed,
% nlat x 240), mm (240x1 calendar-month index), lat.
%
% Haipeng Zhao
% 2026-07-26

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

mm = repmat([1 : 12], 1, 20)';

mldLat = squeeze( nanmean( mld ) );
save( fullfile(OUT_DIR, 'mldLat_monthly_d125'), 'mldLat', 'mm', 'lat' )
fprintf('Saved: %s\n', fullfile(OUT_DIR, 'mldLat_monthly_d125.mat'));
