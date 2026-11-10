% Used for: regrid N2 to 1 deg (Supplementary Table 1)
% Original: Manuscript\MLD\manuscript\nature\NCC\R1\code for revision\correlation_bvf_npp\prj_grid_bvf200_weighted_ot.m (data paths updated to I:\Prj; computation unchanged)
% prj_grid_bvf200_weighted_ot.m
%
% Stage 1 of 3: regrid the new depth-weighted 0-200m BVF (0.083 deg,
% 1441x4320, -60..60 lat only -- I:\Prj\N2_1M_NEMO\bvf_200_depth_weighted\
% <year>\bvf_200_<year>_<month>.nc, field bvf_200_weighted, produced
% 2026-07-27/28 with the np.trapezoid-vs-numpy1.24 fix) to the shared 1x1
% degree grid, mirroring prj_grid008to1G_N2_ot.m's role in the established
% N2 pipeline.
%
% The existing mldG008to1G_idx.mat mapping does NOT apply here -- it was
% built for the OLD N2 pipeline's native grid (4320x2041, confirmed by
% gridIdx numel = 4320*2041), a different domain than this data's
% 1441x4320 (-60..60 lat only) grid. This script builds a fresh
% bin-average mapping instead, using the actual lat/lon coordinates of
% the depth-weighted BVF files directly (simple nearest-1-degree-bin
% mean), matching the same coarse-graining logic (average all fine cells
% falling in each 1deg target cell) as the original.
%
% Target grid: same as the shared 1x1deg grid (lat -89.5..89.5,
% lon -179.5..179.5 in native -180..180 form here since m_map reordering
% only matters for plotting, not for this data stage), taken directly
% from the depth-weighted BVF's own coordinate convention. Cells poleward
% of 60 deg are NaN (source data doesn't cover them).
%
% Output: I:\Prj\N2_1M_NEMO\bvf_200_weighted_one\bvf200_<year>.mat --
% bvf200G1deg (180,360,12) = [lat,lon,month], lat, lon (matches the
% N2G1deg / lat / lon convention of the old pipeline's regridded output).
%
% Haipeng Zhao
% 2026-07-29

clear; clc;
BVF_DIR = 'I:\Prj\N2_1M_NEMO\bvf_200_depth_weighted';
OUT_DIR = 'I:\Prj\N2_1M_NEMO\bvf_200_weighted_one';
if ~exist(OUT_DIR,'dir'), mkdir(OUT_DIR); end

YEARS = 2003:2022;

% Target 1deg grid (matches the shared grid's lat/lon convention exactly)
lat1 = (-89.5:1:89.5)';
lon1 = (-179.5:1:179.5)';
NLAT1 = numel(lat1); NLON1 = numel(lon1);

% Read the fine grid coordinates once (identical across all files)
ref = ncinfo(fullfile(BVF_DIR, '2003', 'bvf_200_2003_01.nc'));
fine_lat = ncread(fullfile(BVF_DIR,'2003','bvf_200_2003_01.nc'), 'lat');
fine_lon = ncread(fullfile(BVF_DIR,'2003','bvf_200_2003_01.nc'), 'lon');
NLATF = numel(fine_lat); NLONF = numel(fine_lon);

% Precompute target-bin index for every fine-grid row/column (1deg bins,
% edges at half-integers matching lat1/lon1 centers)
lat_bin = discretize(fine_lat, -90:1:90);   % 1..180, NaN outside
lon_bin = discretize(fine_lon, -180:1:180); % 1..360

fprintf('Fine grid: %d x %d.  Building fast bin-average indices...\n', NLATF, NLONF);
[LON_BIN, LAT_BIN] = meshgrid(lon_bin, lat_bin);   % both (NLATF,NLONF), matches (lat,lon) field orientation
valid_bin = isfinite(LAT_BIN) & isfinite(LON_BIN);
lin_target = sub2ind([NLAT1, NLON1], LAT_BIN(valid_bin), LON_BIN(valid_bin));

for y = 1:numel(YEARS)
    yyyy = YEARS(y);
    bvf200G1deg = nan(NLAT1, NLON1, 12, 'single');
    for m = 1:12
        f = fullfile(BVF_DIR, num2str(yyyy), sprintf('bvf_200_%d_%02d.nc', yyyy, m));
        v = ncread(f, 'bvf_200_weighted');   % actual on-disk dims (lon,lat) = (4320,1441) per ncdisp -- NOT (lat,lon) as originally assumed here
        v = v';                              % transpose to (lat,lon) = (1441,4320) to match fine_lat/fine_lon convention used to build valid_bin/lin_target
        v = v(valid_bin);
        sums = accumarray(lin_target, v, [NLAT1*NLON1, 1], @(x) mean(x,'omitnan'), NaN);
        bvf200G1deg(:,:,m) = single(reshape(sums, NLAT1, NLON1));
    end
    save(fullfile(OUT_DIR, sprintf('bvf200_%d.mat', yyyy)), 'bvf200G1deg', 'lat1', 'lon1', '-v7.3');
    fprintf('%d done\n', yyyy);
end
fprintf('All years regridded to 1deg: %s\n', OUT_DIR);
