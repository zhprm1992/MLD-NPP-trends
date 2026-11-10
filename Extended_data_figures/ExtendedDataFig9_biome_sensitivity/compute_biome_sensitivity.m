% Used for: biome sensitivity, senMLD_biome_nemo_sig005.mat (Extended Data Fig. 9)
% Original: Manuscript\MLD\manuscript\nature\NCC\R1\code for revision\mld_npp_relationship_biome\compute_biome_sensitivity.m (data paths updated to I:\Prj; computation unchanged)
% compute_biome_sensitivity.m
%
% Stage 2+3 of the biome-scale sensitivity analysis: biome-scale analog
% of Code\npp_cafe\senMLD2_ot.m, replacing "lat row" with "grid cell
% within biome" as the unit that supplies multiple (NPP_trend, MLD_trend)
% pairs per group (confirmed with user 2026-07-30 -- mirrors the lat-row
% logic exactly, just at grid-cell resolution).
%
% For each biome and each of the 4 hemisphere-aware seasons: pool all
% (grid cell x month-in-season) pairs from compute_grid_monthly_trends.m's
% output, keep only pairs where BOTH the NPP and MLD monthly trend have
% p < 0.1 (same threshold as senMLD2_ot.m), then sensitivity =
% median(NPP_trend_pct / MLD_trend_pct) over those pairs.
%
% All 21 biomes are used (no ICE/marginal-sea exclusion) -- consistent
% with Fig. 4's own lat-band version, which pools every latitude row with
% no biome-type filtering. (This differs from the *separate* 14-biome
% grid-cell comparison analysis in ..\mld_npp_relationship_grid, which
% does exclude ICE biomes for a different reason/purpose.)
%
% Biome idx convention: linear index into a (360,180) = [lon,lat] shaped
% array (verified convention, see regress_biome_mld_npp_trends.m header) --
% grid_monthly_trends_mld_npp_nemo.mat's arrays are (180,360,12) =
% [lat,lon,month], so each (lat,lon) slice must be transposed to (lon,lat)
% before applying idx.
%
% Output: senMLD_biome_nemo.mat (this folder) -- senMLD (nbiome,4,3) =
% [median, std, n] per biome/season, senPt (nbiome x 4 cell array) of raw
% [npp_trend, mld_trend] pairs (mld set to 0 for non-significant pairs,
% matching senMLD2_ot.m's convention for the scatter-overlay plot),
% biome_names, season_names.
%
% Haipeng Zhao
% 2026-07-30

function compute_biome_sensitivity(sig_level, out_suffix)
% compute_biome_sensitivity(sig_level, out_suffix)
%   sig_level  : significance threshold for both NPP and MLD trend p-values
%                (e.g. 0.1 or 0.05). Default 0.1.
%   out_suffix : appended to the output filename, e.g. 'sig010', 'sig005'.
%                Default derived from sig_level.
%
% Excludes all Arctic biomes (NP ICE + the 5 Arctic marginal seas) per
% user request 2026-07-30 -- SO ICE (Antarctic) is kept.

if nargin < 1 || isempty(sig_level), sig_level = 0.1; end
if nargin < 2 || isempty(out_suffix), out_suffix = sprintf('sig%03d', round(sig_level*100)); end

OUT_DIR = fileparts(mfilename('fullpath'));
d = load(fullfile(OUT_DIR, 'grid_monthly_trends_mld_npp_nemo.mat'));

mmN = [3 4 5; 6 7 8; 9 10 11; 12 1 2];   % [Spring, Summer, Autumn, Winter], lat >= 0
mmS = [9 10 11; 12 1 2; 3 4 5; 6 7 8];   % lat < 0
season_names = {'Spring','Summer','Autumn','Winter'};

% All polar/ICE-regime biomes excluded (NP ICE, SO ICE) + Arctic marginal
% seas (Baffin Bay, Greenland Sea, Hudson Bay, Kara Sea, Laptev Sea), per
% user request 2026-07-30 (no polar regime considered).
biomes = {'A EQU','EP EQU','WP EQU','NA STPS','SA STPS','NP STPS','SP STPS','IND STPS', ...
    'NA STSS','NP STSS','SO STSS','NA SPSS','NP SPSS','SO SPSS'};
NBIOME = numel(biomes);
BIOME_DIR = 'I:\Prj\MLD_1M_NEMO\biomes_1M_1dg\mld\2003';

lat = d.lat(:); NLAT = numel(lat); NLON = numel(d.lon);
isN = lat >= 0;  % (180,1), broadcasts across lon after transpose to (lon,lat)

senMLD = nan(NBIOME, 4, 3);
senPt = cell(NBIOME, 4);
biome_n_pixels = nan(NBIOME,1);

for bi = 1:NBIOME
    bname = biomes{bi};
    bd = load(fullfile(BIOME_DIR, sprintf('mld_%s.mat', bname)), 'idx');
    idx = round(bd.idx(:));
    biome_n_pixels(bi) = numel(idx);
    fprintf('%-16s n_pixels=%d\n', bname, numel(idx));

    for s = 1:4
        npp_pts = []; mld_pts = []; npp_sig_pts = []; mld_sig_pts = [];
        months_this_season_N = mmN(s,:);
        months_this_season_S = mmS(s,:);

        for mi = 1:3
            % NPP/MLD trend+p slices for this calendar month, hemisphere-aware:
            % northern latitudes use months_this_season_N(mi), southern use
            % months_this_season_S(mi) -- both are (180,360) = [lat,lon].
            npp_trend_N = d.npp_trend_pct(:,:,months_this_season_N(mi));
            npp_p_N     = d.npp_trend_p(:,:,months_this_season_N(mi));
            npp_trend_S = d.npp_trend_pct(:,:,months_this_season_S(mi));
            npp_p_S     = d.npp_trend_p(:,:,months_this_season_S(mi));
            mld_trend_N = d.mld_trend_pct(:,:,months_this_season_N(mi));
            mld_p_N     = d.mld_trend_p(:,:,months_this_season_N(mi));
            mld_trend_S = d.mld_trend_pct(:,:,months_this_season_S(mi));
            mld_p_S     = d.mld_trend_p(:,:,months_this_season_S(mi));

            npp_trend = nan(NLAT,NLON); npp_trend(isN,:) = npp_trend_N(isN,:); npp_trend(~isN,:) = npp_trend_S(~isN,:);
            npp_p     = nan(NLAT,NLON); npp_p(isN,:)     = npp_p_N(isN,:);     npp_p(~isN,:)     = npp_p_S(~isN,:);
            mld_trend = nan(NLAT,NLON); mld_trend(isN,:) = mld_trend_N(isN,:); mld_trend(~isN,:) = mld_trend_S(~isN,:);
            mld_p     = nan(NLAT,NLON); mld_p(isN,:)     = mld_p_N(isN,:);     mld_p(~isN,:)     = mld_p_S(~isN,:);

            % transpose (lat,lon) -> (lon,lat) to match idx convention
            npp_trend_t = npp_trend'; npp_p_t = npp_p';
            mld_trend_t = mld_trend'; mld_p_t = mld_p';

            npp_pts = [npp_pts; npp_trend_t(idx)];
            mld_pts = [mld_pts; mld_trend_t(idx)];
            npp_sig_pts = [npp_sig_pts; npp_p_t(idx)];
            mld_sig_pts = [mld_sig_pts; mld_p_t(idx)];
        end

        valid = ~isnan(npp_pts);
        npp_pts = npp_pts(valid); mld_pts = mld_pts(valid);
        npp_sig_pts = npp_sig_pts(valid); mld_sig_pts = mld_sig_pts(valid);

        idx_sig = (npp_sig_pts < sig_level) & (mld_sig_pts < sig_level);
        sen = npp_pts(idx_sig) ./ mld_pts(idx_sig);
        sen = sen(~isinf(sen) & ~isnan(sen));
        if ~isempty(sen)
            senMLD(bi, s, 1) = nanmedian(sen);
            senMLD(bi, s, 2) = nanstd(sen);
            senMLD(bi, s, 3) = numel(sen);
        end

        mld_pts_plot = mld_pts;
        mld_pts_plot(~idx_sig) = 0;
        senPt{bi, s} = [npp_pts, mld_pts_plot];
    end
end

out_file = fullfile(OUT_DIR, sprintf('senMLD_biome_nemo_%s.mat', out_suffix));
save(out_file, 'senMLD','senPt','biomes','biome_n_pixels','season_names','sig_level');
fprintf('Saved: %s\n', out_file);
end
