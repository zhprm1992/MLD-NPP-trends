% Used for: per-grid-cell monthly trends of VGPM and CAFE NPP (Fig. R4)
% Original: Manuscript\MLD\manuscript\nature\NCC\R2\code for revision\interannual_biome\grid_monthly_trends_vgpm_cafe.m (data paths updated to I:\Prj; computation unchanged)
% grid_monthly_trends_vgpm_cafe.m
% Per 1x1 grid cell and calendar month, 2003-2022 NPP trend (% and Mann-Kendall p) for VGPM and CAFE-NEMO.
% Exact copy of the logic in R1\code for revision\mld_npp_relationship_biome\compute_grid_monthly_trends.m
% (ktaub(datain, 0.9999, 1), pct = (fit_end - fit_start) / fit_start * 100, all 20 years finite), NPP only:
% MLD trends are reused from the R1 file. Output: grid_monthly_trends_npp_vgpm_cafe.mat (this folder).
% READ ONLY on inputs. Haipeng Zhao / Claude, 2026-11-05 (R2)

clear; clc;
OUT = fileparts(mfilename('fullpath'));
diary(fullfile(OUT, 'grid_monthly_trends_vgpm_cafe_log.txt')); diary on
ROOT = 'D:\Users\zhp19\OneDrive - Indiana State University\Research';
if ~isfolder(ROOT), ROOT = 'E:\OneDrive - Indiana State University\Research'; end
addpath(fullfile(ROOT, 'MATLAB toolbox', 'ktaub'));
YEARS = 2003:2022; NYR = numel(YEARS); x = (1:NYR)';
src = {'VGPM', 'I:\Prj\npp_vgpm_1M_1dg\grid_npp_vgpm_1M_1dg_mat\npp_vgpm_YYYY.mat';
       'CAFE', 'I:\Prj\npp_cafe_1M_1dg\grid_npp_nemo_1M_1dg_mat\npp_YYYY.mat'};
for k = 1:2
    fprintf('%s: loading ... %s\n', src{k,1}, datestr(now));
    X = nan(180, 360, 12, NYR, 'single');
    for y = 1:NYR, d = load(strrep(src{k,2}, 'YYYY', num2str(YEARS(y))), 'nppG1dg'); X(:,:,:,y) = d.nppG1dg; end
    pct = nan(180, 360, 12, 'single'); pv = pct;
    for m = 1:12
        tm = tic;
        for i = 1:180
            for j = 1:360
                v = double(squeeze(X(i,j,m,:)));
                if any(isnan(v)), continue; end
                [~,~,~,sig,~,~,~,~,~,pos,~] = ktaub([x v], 0.9999, 1);
                if ~isnan(pos)
                    if pos(1,2) == 0, pct(i,j,m) = 0; else, pct(i,j,m) = single((pos(2,2) - pos(1,2)) / pos(1,2) * 100); end
                    pv(i,j,m) = sig;
                end
            end
        end
        fprintf('  %s month %d done in %.0fs\n', src{k,1}, m, toc(tm));
    end
    T.([lower(src{k,1}) '_trend_pct']) = pct; T.([lower(src{k,1}) '_trend_p']) = pv;
end
save(fullfile(OUT, 'grid_monthly_trends_npp_vgpm_cafe.mat'), '-struct', 'T', '-v7.3');
fprintf('Saved %s\n', datestr(now)); diary off
