% Used for: GLMM per model and season with latitude band as group (Fig. R2)
% Original: Manuscript\MLD\manuscript\nature\NCC\R2\code for revision\interannual_fig4\glmm_band_group.m (data paths updated to I:\Prj; computation unchanged)
% glmm_band_group.m
% Reviewer 2 comment 1, GLMM with latitudinal band as the group (HZ decision 2026-11-05, see revision_log_R2.md).
% Data: detrended ln NPP and ln MLD-NEMO of the Fig. 4 samples (lat row x month with paired significant trends,
% p < 0.1), saved by interannual_fig4_samples.m. One model per NPP model x season, all six bands together:
%   main : lnNPP ~ lnMLD + (1 + lnMLD | band)
%   check: lnNPP ~ lnMLD + (1 + lnMLD | band) + (lnMLD - 1 | sample)      (samples nested in bands)
% Band slope = fixed slope + band random effect; approximate 95% CI from sqrt(SE_fixed^2 + SEPred^2).
% READ ONLY on inputs. Haipeng Zhao / Claude, 2026-11-05 (R2)

clear; clc;
OUT = fileparts(mfilename('fullpath'));
sfx = '_pct';                                 % '' = log version, '_pct' = percent anomaly version
diary(fullfile(OUT, ['glmm_band_group' sfx '_log.txt'])); diary on
fprintf('Start %s\n', datestr(now));
load(fullfile(OUT, ['interannual_fig4_samples' sfx '.mat']), 'P', 'S');
mmN = [3 4 5; 6 7 8; 9 10 11; 12 1 2];  mmS = [9 10 11; 12 1 2; 3 4 5; 6 7 8];
seasName = {'Spring','Summer','Autumn','Winter'};
latBand = [60 40; 40 20; 20 0; 0 -20; -20 -40; -40 -60];
bandName = {'40-60N','20-40N','0-20N','0-20S','20-40S','40-60S'};
models = {'CbPM','VGPM','CAFE'};
R = {};
for k = 1:3
    ok = S.(models{k}).ok; latF = S.(models{k}).latF; ser = P.(models{k}).series;
    for s = 1:4
        rows = []; nS = zeros(6,1); sid = 0;
        for u = 1:6
            [r, c] = sampleIdx(ok, latF, latBand(u,:), mmN(s,:), mmS(s,:));
            for q = 1:numel(r)
                z = ser{r(q), c(q)};
                if isempty(z), continue; end
                sid = sid + 1; nS(u) = nS(u) + 1;
                rows = [rows; z repmat([u sid], size(z,1), 1)]; %#ok<AGROW>
            end
        end
        T = table(rows(:,1), rows(:,2), categorical(rows(:,3)), categorical(rows(:,4)), 'VariableNames', {'lnNPP','lnMLD','band','sample'});
        for v = 1:2
            if v == 1, f = 'lnNPP ~ lnMLD + (1 + lnMLD | band)'; tag = 'main';
            else, f = 'lnNPP ~ lnMLD + (1 + lnMLD | band) + (lnMLD - 1 | sample)'; tag = 'nested'; end
            try
                m = fitlme(T, f); ci = coefCI(m); b0 = m.Coefficients.Estimate(2); se0 = m.Coefficients.SE(2);
                [~, ~, st] = randomEffects(m);
                isB = strcmp(st.Group, 'band') & strcmp(st.Name, 'lnMLD');
                lev = str2double(st.Level(isB)); bb = st.Estimate(isB); sp = st.SEPred(isB);
                R(end+1,:) = {models{k}, seasName{s}, tag, 'All bands', sum(nS), b0, ci(2,1), ci(2,2)}; %#ok<SAGROW>
                for u = 1:6
                    j = find(lev == u);
                    if nS(u) == 0 || isempty(j), est = NaN; lo = NaN; hi = NaN;
                    else, est = b0 + bb(j); h = 1.96 * sqrt(se0^2 + sp(j)^2); lo = est - h; hi = est + h; end
                    R(end+1,:) = {models{k}, seasName{s}, tag, bandName{u}, nS(u), est, lo, hi}; %#ok<SAGROW>
                end
            catch e
                fprintf('  LME fail %s %s %s: %s\n', models{k}, seasName{s}, tag, e.message);
            end
        end
    end
end
Tab = cell2table(R, 'VariableNames', {'npp','season','model','band','n_samples','slope','ci_lo','ci_hi'});
writetable(Tab, fullfile(OUT, ['glmm_band_group' sfx '.csv']));
save(fullfile(OUT, ['glmm_band_group' sfx '.mat']), 'Tab');

for tag = {'main','nested'}
    fprintf('\n===== %s model: slope [95%% CI] (n samples); rows = band, cols Sp Su Au Wi =====\n', tag{1});
    for k = 1:3
        fprintf('-- %s\n', models{k});
        for un = [{'All bands'}, bandName]
            fprintf('  %-10s', un{1});
            for s = 1:4
                q = Tab(strcmp(Tab.npp, models{k}) & strcmp(Tab.model, tag{1}) & strcmp(Tab.band, un{1}) & strcmp(Tab.season, seasName{s}), :);
                if isempty(q), fprintf('   %-24s', 'fail'); continue; end
                fl = ' '; if q.ci_hi < 0, fl = '-'; elseif q.ci_lo > 0, fl = '+'; end
                fprintf('  %5.2f[%5.2f,%5.2f]%s(%3d)', q.slope, q.ci_lo, q.ci_hi, fl, q.n_samples);
            end
            fprintf('\n');
        end
    end
    fprintf('Band-seasons (24): CI<0 / CI>0 / spans 0 / no samples\n');
    for k = 1:3
        t = Tab(strcmp(Tab.npp, models{k}) & strcmp(Tab.model, tag{1}) & ~strcmp(Tab.band, 'All bands'), :);
        fprintf('  %s: %2d / %2d / %2d / %2d\n', models{k}, sum(t.ci_hi < 0), sum(t.ci_lo > 0), sum(t.ci_lo <= 0 & t.ci_hi >= 0), sum(t.n_samples == 0));
    end
end
fprintf('\nEnd %s\n', datestr(now)); diary off

function [r, c] = sampleIdx(ok, latF, band, mN, mS)
    r = []; c = [];
    inB = latF < band(1) & latF > band(2);
    for h = 1:2
        if h == 1, rows = inB & latF > 0; mm = mN; else, rows = inB & latF < 0; mm = mS; end
        [rr, cc] = find(ok(:, mm) & rows); r = [r; rr]; c = [c; mm(cc)']; %#ok<AGROW>
    end
end
