% Used for: ensemble-mean NPP (VGPM, CAFE, CbPM) per grid cell (SI Figs. S3-S5)
% Original: Code\npp_ensemble\grid_npp_ensemble_1M_1dg_allProducts.m (data paths updated to I:\Prj; computation unchanged)
% Output ensemble-mean NPP from VGPM, CAFE, and CbPM, for all 3 MLD
% products (NEMO, ECCO, HYCOM) -- generalizes grid_npp_ensemble_1M_1dg.m
% (NEMO-only) to loop over products, and writes output in the same
% folder/file convention as CAFE (npp_cafe_1M_1dg\grid_npp_<product>_
% 1M_1dg_mat\npp_<year>.mat, variable nppG1dg) -- i.e. subfolder names
% keep the "_1M_1dg_mat" suffix, matching CAFE rather than CbPM's plain
% grid_npp_<product> (no suffix). File naming (npp_<year>.mat) matches
% both conventions, which are identical on that point.
%
% 2026-07-15: renamed grid_npp_<product> -> grid_npp_<product>_1M_1dg_mat
% after first generating with the CbPM-style (no-suffix) names; old
% NEMO-only grid_npp_ensemble_1M_1dg_mat folder deleted as a duplicate
% (verified byte-identical to the new grid_npp_nemo_1M_1dg_mat first).
%
% Ensemble = pixel-wise, month-wise mean across the three algorithms
% (omitting NaN where an individual model has no data), same method as
% the NEMO-only version.
%
% VGPM does not depend on MLD product (single universal grid, no per-
% product forcing) -- the same VGPM file is reused for all 3 ensembles.
% CAFE and CbPM each have NEMO/ECCO/HYCOM-forced variants and are paired
% product-for-product.
%
% NOTE (carried over from the NEMO-only version): eVGPM deliberately
% excluded (shares Pb_opt with VGPM, would double-count that bias). CAFE
% and CbPM currently only have 2003-2022 (20 years) vs. VGPM's
% 2003-2023 -- ensemble covers the overlapping 2003-2022 only.
%
% Haipeng Zhao
% 2026-07-15

clear
clc

products = {'nemo','ecco','hycom'};

PATH_IN_vgpm = 'I:\Prj\npp_vgpm_1M_1dg\grid_npp_vgpm_1M_1dg_mat';   % same for all products
PATH_IN_chl  = 'I:\Prj\npp_cafe_1M_1dg\grid_chl_1M_1dg_mat';        % for lat/lon only

PATH_OT_base = 'I:\Prj\npp_ensemble_1M_1dg';

file_chl = dir( [PATH_IN_chl, '\', '*.mat'] );
load( [file_chl(1).folder, '\', file_chl(1).name], 'lat', 'lon' );

for pi = 1:numel(products)
    prod = products{pi};
    fprintf('=== Product: %s ===\n', upper(prod));

    PATH_IN_cafe = fullfile('I:\Prj\npp_cafe_1M_1dg', ['grid_npp_' prod '_1M_1dg_mat']);
    PATH_IN_cbpm = fullfile('I:\Prj\npp_cbpm_1M_1dg_zno3mMU', ['grid_npp_' prod]);
    PATH_OT_npp  = fullfile(PATH_OT_base, ['grid_npp_' prod '_1M_1dg_mat']);
    if ~exist(PATH_OT_npp, 'dir'), mkdir(PATH_OT_npp); end

    for yyyy = 2003 : 2022
        yyyyStr = num2str( yyyy );

        d_vgpm = load( fullfile( PATH_IN_vgpm, ['npp_vgpm_', yyyyStr, '.mat'] ), 'nppG1dg' );
        d_cafe = load( fullfile( PATH_IN_cafe, ['npp_', yyyyStr, '.mat'] ), 'nppG1dg' );
        d_cbpm = load( fullfile( PATH_IN_cbpm, ['npp_', yyyyStr, '.mat'] ), 'nppG1dg' );

        npp_stack = cat( 4, d_vgpm.nppG1dg, d_cafe.nppG1dg, d_cbpm.nppG1dg ); % [180, 360, 12, 3]
        nppG1dg = mean( npp_stack, 4, 'omitnan' );
        n_models = sum( ~isnan( npp_stack ), 4 ); % QC: number of models contributing per pixel/month

        filename_ot = ['npp_', yyyyStr];   % matches CbPM's npp_<year>.mat convention (no '_ensemble' infix)
        save( fullfile( PATH_OT_npp, filename_ot ), 'nppG1dg', 'n_models', 'lat', 'lon', '-v7.3' )
        fprintf('  saved %s\n', fullfile(PATH_OT_npp, filename_ot));
    end
end

fprintf('\nDone. Output layout mirrors CAFE (npp_cafe_1M_1dg\\grid_npp_<product>_1M_1dg_mat):\n');
for pi = 1:numel(products)
    fprintf('  %s\n', fullfile(PATH_OT_base, ['grid_npp_' products{pi} '_1M_1dg_mat'], 'npp_<year>.mat'));
end
