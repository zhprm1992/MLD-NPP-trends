% Used for: ensemble biome sensitivity, senMLD_biome_ensemble_sig005.mat (Supplementary Fig. 5)
% Original: Manuscript\MLD\manuscript\nature\NCC\R1\code for revision\mld_npp_relationship_biome\compute_biome_sensitivity_ensemble.m (data paths updated to I:\Prj; computation unchanged)
% compute_biome_sensitivity_ensemble.m
%
% Ensemble-NPP counterpart to compute_biome_sensitivity.m. Identical
% method (pool grid-cell x month-in-season (NPP_trend, MLD_trend) pairs
% per biome/season, keep only pairs where BOTH trends have p < sig_level,
% sensitivity = median(NPP_trend_pct / MLD_trend_pct) over those pairs),
% reading grid_monthly_trends_mld_npp_ENSEMBLE_nemo.mat (ensemble NPP)
% instead of the CbPM file. Same 14-biome list (excludes ICE + Arctic
% marginal seas), same BIOME_DIR, same idx convention.
%
% Output: senMLD_biome_ensemble_<suffix>.mat (this folder) -- same
% structure as senMLD_biome_nemo_<suffix>.mat.
%
% Haipeng Zhao
% 2026-09-12

function compute_biome_sensitivity_ensemble(sig_level, out_suffix)
if nargin < 1 || isempty(sig_level), sig_level = 0.1; end
if nargin < 2 || isempty(out_suffix), out_suffix = sprintf('sig%03d', round(sig_level*100)); end

OUT_DIR = fileparts(mfilename('fullpath'));
d = load(fullfile(OUT_DIR, 'grid_monthly_trends_mld_npp_ensemble_nemo.mat'));

mmN = [3 4 5; 6 7 8; 9 10 11; 12 1 2];   % [Spring, Summer, Autumn, Winter], lat >= 0
mmS = [9 10 11; 12 1 2; 3 4 5; 6 7 8];   % lat < 0
season_names = {'Spring','Summer','Autumn','Winter'};

biomes = {'A EQU','EP EQU','WP EQU','NA STPS','SA STPS','NP STPS','SP STPS','IND STPS', ...
    'NA STSS','NP STSS','SO STSS','NA SPSS','NP SPSS','SO SPSS'};
NBIOME = numel(biomes);
BIOME_DIR = 'I:\Prj\MLD_1M_NEMO\biomes_1M_1dg\mld\2003';

lat = d.lat(:); NLAT = numel(lat); NLON = numel(d.lon);
isN = lat >= 0;

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

out_file = fullfile(OUT_DIR, sprintf('senMLD_biome_ensemble_%s.mat', out_suffix));
save(out_file, 'senMLD','senPt','biomes','biome_n_pixels','season_names','sig_level');
fprintf('Saved: %s\n', out_file);
end
