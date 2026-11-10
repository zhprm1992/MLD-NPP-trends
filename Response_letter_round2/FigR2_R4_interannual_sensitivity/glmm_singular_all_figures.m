% Used for: flags singular GLMM fits (no significance mark, Figs. R2, R4)
% Original: Manuscript\MLD\manuscript\nature\NCC\R2\code for revision\interannual_fig4\glmm_singular_all_figures.m (data paths updated to I:\Prj; computation unchanged)
% glmm_singular_all_figures.m  (2026-11-08, same rule as glmm_singular_check.m, applied to every figure that plots group slopes)
% Singular = group slopes identical within a fit (range < 1e-3 among groups with samples, >= 2 groups).
% Figures 01/03 (glmm_band_group_pct), 04 (lag_check_fig4, lag0 and season rows), 02 (interannual_biome_pct). READ ONLY.
clear; clc; OUT = fileparts(mfilename('fullpath')); R2 = fileparts(OUT);
diary(fullfile(OUT, 'glmm_singular_all_figures_log.txt')); diary on
seasName = {'Spring','Summer','Autumn','Winter'}; models = {'CbPM','VGPM','CAFE'}; SING = {};
G = load(fullfile(OUT, 'glmm_band_group_pct.mat')); A = G.Tab(strcmp(G.Tab.model, 'main') & ~strcmp(G.Tab.band, 'All bands') & G.Tab.n_samples > 0, :);
SING = [SING; chk(A, 'slope', 'fig01/03', models, seasName, {}, '')];
L = load(fullfile(OUT, 'lag_check_fig4.mat'), 'Tab'); disp(L.Tab.Properties.VariableNames)
for lg = {'lag0', 'season'}
    A = L.Tab(strcmp(L.Tab.lag, lg{1}) & isfinite(L.Tab.glmm), :);
    SING = [SING; chk(A, 'glmm', ['fig04 ' lg{1}], models, seasName, {}, lg{1})];
end
B = load(fullfile(R2, 'interannual_biome', 'interannual_biome_pct.mat'), 'Tab'); disp(B.Tab.Properties.VariableNames)
SING = [SING; chk(B.Tab, '', 'fig02', models, seasName, {}, '')];
S = cell2table(SING, 'VariableNames', {'figure','npp','season','lag','n_groups','range','singular'});
save(fullfile(OUT, 'glmm_singular_all_figures.mat'), 'S'); writetable(S, fullfile(OUT, 'glmm_singular_all_figures.csv'));
fprintf('\nSINGULAR fits:\n'); disp(S(S.singular, :)); diary off
function R = chk(A, col, fig, models, seasName, ~, lg)
    R = {}; vn = A.Properties.VariableNames;
    if isempty(col), col = vn{find(ismember(vn, {'slope','glmm','band_slope','biome_slope'}), 1)}; end
    nc = vn{find(ismember(vn, {'npp','model'}), 1)};
    gc = vn{find(ismember(vn, {'band','biome'}), 1)};
    for k = 1:3
        for s = 1:4
            w = strcmp(A.(nc), models{k}) & strcmp(A.season, seasName{s}) & isfinite(A.(col));
            if ismember('n_samples', vn), w = w & A.n_samples > 0; end
            if ismember(gc, vn) && iscell(A.(gc)), w = w & ~ismember(A.(gc), {'All bands','All biomes','All'}); end
            v = A.(col)(w); if numel(v) < 2, continue; end
            rg = max(v) - min(v); sg = rg < 1e-3;
            fprintf('%-12s %-5s %-7s groups %2d  slope range %.2e %s\n', fig, models{k}, seasName{s}, numel(v), rg, repmat('SINGULAR', 1, sg));
            R(end+1,:) = {fig, models{k}, seasName{s}, lg, numel(v), rg, sg}; %#ok<AGROW>
        end
    end
end
