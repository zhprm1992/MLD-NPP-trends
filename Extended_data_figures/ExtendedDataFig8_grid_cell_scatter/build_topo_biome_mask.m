% Used for: open-ocean mask topo_biome_mask_nemo.mat (Extended Data Figs. 4, 7 and 8, Supplementary Fig. 4)
% Original: Manuscript\MLD\manuscript\nature\NCC\R1\code for revision\mld_mechanism\build_topo_biome_mask.m (data paths updated to I:\Prj; computation unchanged)
% build_topo_biome_mask.m
%
% Reusable topo+biome validity mask, replicating the EXACT masking used
% by the real Figure 4 pipeline's lat-row products (Code\mld_NEMO\
% monthlyLat_ot.m, Code\npp_cbpm\monthlyLat_nemo_ot.m): a grid cell
% counts only if (a) it belongs to one of the 21 established biomes and
% (b) topo <= -1000m (deep ocean only, shallow/coastal cells excluded).
%
% Found 2026-08-24: my own lat-row scripts this session (compute_latrow_
% monthly_trends_{npp,parmld_mld,igfun_mld}.m, compute_grid_light_
% attribution_mld.m's lat-row companion, build_nitrate_latmonthly_
% climatology.m) all used a NAIVE unmasked longitude mean instead of
% this established convention -- confirmed by a direct mismatch against
% the real nppLat_monthly_nemo.mat (1039 vs 1108 at a spot-checked
% cell). This script builds the correct mask once, for all of those to
% be rebuilt against.
%
% idx convention: linear index into a (360,180) = [lon,lat] shaped array
% (established convention, verified repeatedly this project).
%
% Output: topo_biome_mask_nemo.mat (this folder) -- valid_lonlat (360x180
% logical, true = counts), lat (180,1), lon (360,1)
%
% SUPERSEDED 2026-08-24: this array-mask approach has an unresolved bug
% (root-caused, not fixed here) -- applying valid_lonlat to a raw grid
% variable overcounts relative to the real biome-split products, because
% the true biomes_<var>_nemo files carry extra per-biome QC NaN-masking
% beyond simple topo+idx membership (confirmed: at lat row 21 Jan 2003,
% 128 non-NaN cells in the real biome-split file vs 150 at the same
% idx/topo-passed cells pulled from the raw grid). Use
% biome_split_zonal_mean.m (reads biome-split files directly, exact
% match validated) or latrow_topo_biome_mean.m (raw-grid fallback for
% variables with no biome-split product, documented approximation)
% instead -- both in this folder.
%
% Haipeng Zhao
% 2026-08-24

clear; clc;
addpath('D:\Users\zhp19\OneDrive - Indiana State University\Research\MATLAB toolbox\m_map');  % m_topo.mat

BIOME_DIR = 'I:\Prj\MLD_1M_NEMO\biomes_1M_1dg\mld\2003';
datafile = dir(fullfile(BIOME_DIR, '*.mat'));
fprintf('Found %d biome files\n', numel(datafile));

ref = load(fullfile(BIOME_DIR, datafile(1).name), 'lat','lon');
lon = ref.lon; lat = ref.lat;
u1 = length(lon); u2 = length(lat);

d = load('m_topo.mat');
topo = d.topo';

valid_lonlat = false(u1, u2);
for j = 1:numel(datafile)
    bd = load(fullfile(BIOME_DIR, datafile(j).name), 'idx');
    for i = 1:length(bd.idx)
        [I, J] = ind2sub([u1, u2], bd.idx(i));
        if topo(I, J) > -1000
            continue
        end
        valid_lonlat(I, J) = true;
    end
end

fprintf('Valid (deep-ocean, biome-member) cells: %d of %d (%.1f%%)\n', ...
    sum(valid_lonlat(:)), numel(valid_lonlat), 100*sum(valid_lonlat(:))/numel(valid_lonlat));

OUT_DIR = fileparts(mfilename('fullpath'));
save(fullfile(OUT_DIR, 'topo_biome_mask_nemo.mat'), 'valid_lonlat', 'lat', 'lon');
fprintf('Saved: %s\n', fullfile(OUT_DIR, 'topo_biome_mask_nemo.mat'));
