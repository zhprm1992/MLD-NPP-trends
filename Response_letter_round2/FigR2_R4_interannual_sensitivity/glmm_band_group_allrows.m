% Used for: interannual GLMM on all latitude rows and months, no significance filter (Fig. R2 bottom row)
% Original: Manuscript\MLD\manuscript\nature\NCC\R2\code for revision\interannual_fig4\glmm_band_group_allrows.m (data paths updated to I:\Prj; computation unchanged)
% glmm_band_group_allrows.m
% Internal check for Nicolas's comment (2026-11-10, see revision_log_R2.md): is the interannual sensitivity of
% Figure R2 representative beyond the Fig. 4 samples? Same series, transform and model as Figure R2
% (interannual_fig4_samples.m with TR = 'pct' and glmm_band_group.m, main model), with ONE change: no
% significance filter. Every 1 deg latitude row in 60S-60N and every calendar month with 20 finite positive
% values is a sample. Model per NPP model x season: lnNPP ~ lnMLD + (1 + lnMLD | band) on detrended percent
% anomalies (variable names kept from the Figure R2 script). READ ONLY on inputs. Haipeng Zhao / Claude, 2026-11-10.

clear; clc;
OUT = fileparts(mfilename('fullpath'));
diary(fullfile(OUT, 'glmm_band_group_allrows_log.txt')); diary on
fprintf('Start %s\n', datestr(now));
tf = @(v) (v - mean(v)) / mean(v) * 100;
mmN = [3 4 5; 6 7 8; 9 10 11; 12 1 2];  mmS = [9 10 11; 12 1 2; 3 4 5; 6 7 8];
seasName = {'Spring','Summer','Autumn','Winter'};
latBand = [60 40; 40 20; 20 0; 0 -20; -20 -40; -40 -60];
bandName = {'40-60N','20-40N','0-20N','0-20S','20-40S','40-60S'};
src = {'CbPM', 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\nppLat_monthly_nemo.mat';
       'VGPM', 'I:\Prj\npp_vgpm_1M_1dg\nppLat_monthly_vgpm.mat';
       'CAFE', 'I:\Prj\npp_cafe_1M_1dg\nppLat_monthly_nemo.mat'};
dm = load('I:\Prj\MLD_1M_NEMO\mldLat_monthly_nemo.mat');
latS = double(dm.lat(:)); keepM = latS < 60 & latS > -60; mldLat = double(dm.mldLat(keepM, :)); latS = latS(keepM);
mmM = dm.mm(:)';
R = {};
for k = 1:3
    dn = load(src{k,2});
    latN = double(dn.lat(:)); keep = latN < 60 & latN > -60; nppLat = double(dn.nppLat(keep, :)); latN = latN(keep);
    if numel(dn.mm) > numel(mmM), nppLat = nppLat(:, 1:numel(mmM)); dn.mm = dn.mm(1:numel(mmM)); end   % VGPM 2003-2023 -> 2003-2022
    assert(isequal(latN, latS), 'lat mismatch'); assert(isequal(dn.mm(:)', mmM), 'month axis mismatch');
    for s = 1:4
        rows = []; nS = zeros(6,1); sid = 0;
        for u = 1:6
            inB = find(latS < latBand(u,1) & latS > latBand(u,2));
            for q = inB'
                if latS(q) > 0, mm = mmN(s,:); else, mm = mmS(s,:); end
                for j = mm
                    y = nppLat(q, mmM == j)'; x = mldLat(q, mmM == j)';
                    if numel(y) ~= 20 || any(~isfinite([y; x])) || any([y; x] <= 0), continue; end
                    sid = sid + 1; nS(u) = nS(u) + 1;
                    rows = [rows; detrend(tf(y)) detrend(tf(x)) repmat([u sid], 20, 1)]; %#ok<AGROW>
                end
            end
        end
        T = table(rows(:,1), rows(:,2), categorical(rows(:,3)), 'VariableNames', {'lnNPP','lnMLD','band'});
        try
            m = fitlme(T, 'lnNPP ~ lnMLD + (1 + lnMLD | band)'); ci = coefCI(m);
            b0 = m.Coefficients.Estimate(2); se0 = m.Coefficients.SE(2);
            [~, ~, st] = randomEffects(m);
            isB = strcmp(st.Group, 'band') & strcmp(st.Name, 'lnMLD');
            lev = str2double(st.Level(isB)); bb = st.Estimate(isB); sp = st.SEPred(isB);
            R(end+1,:) = {src{k,1}, seasName{s}, 'All bands', sum(nS), b0, ci(2,1), ci(2,2)}; %#ok<SAGROW>
            for u = 1:6
                j = find(lev == u);
                if nS(u) == 0 || isempty(j), est = NaN; lo = NaN; hi = NaN;
                else, est = b0 + bb(j); h = 1.96 * sqrt(se0^2 + sp(j)^2); lo = est - h; hi = est + h; end
                R(end+1,:) = {src{k,1}, seasName{s}, bandName{u}, nS(u), est, lo, hi}; %#ok<SAGROW>
            end
        catch e
            fprintf('  LME fail %s %s: %s\n', src{k,1}, seasName{s}, e.message);
        end
    end
end
A = cell2table(R, 'VariableNames', {'npp','season','band','n_samples','slope','ci_lo','ci_hi'});
writetable(A, fullfile(OUT, 'glmm_band_group_allrows.csv'));
save(fullfile(OUT, 'glmm_band_group_allrows.mat'), 'A');

% comparison with Figure R2 (Fig. 4 samples, main model)
F = load(fullfile(OUT, 'glmm_band_group_pct.mat')); F = F.Tab; F = F(strcmp(F.model, 'main'), :);
models = src(:,1)';
fprintf('\nALL ROWS vs FIG. 4 SAMPLES: band slope [95%% CI] (n samples); flag - / + = CI excludes 0\n');
for k = 1:3
    fprintf('-- %s\n', models{k});
    for un = [{'All bands'}, bandName]
        fprintf('  %-9s', un{1});
        for s = 1:4
            a = A(strcmp(A.npp, models{k}) & strcmp(A.band, un{1}) & strcmp(A.season, seasName{s}), :);
            f = F(strcmp(F.npp, models{k}) & strcmp(F.band, un{1}) & strcmp(F.season, seasName{s}), :);
            fa = fl(a); ff = fl(f);
            fprintf(' | %5.2f%s(%3d) vs %5.2f%s(%3d)', a.slope, fa, a.n_samples, f.slope, ff, f.n_samples);
        end
        fprintf('\n');
    end
end
fprintf('\nSummary per model over 24 band-seasons:\n');
for k = 1:3
    a = A(strcmp(A.npp, models{k}) & ~strcmp(A.band, 'All bands'), :);
    f = F(strcmp(F.npp, models{k}) & ~strcmp(F.band, 'All bands'), :);
    [~, ia, ib] = intersect(strcat(a.season, a.band), strcat(f.season, f.band));
    both = isfinite(a.slope(ia)) & isfinite(f.slope(ib));
    fprintf('  %s all rows: CI<0 %2d, CI>0 %2d, spans 0 %2d | Fig.4 samples: CI<0 %2d, CI>0 %2d, spans 0 %2d, none %2d | same sign %d of %d\n', ...
        models{k}, sum(a.ci_hi < 0), sum(a.ci_lo > 0), sum(a.ci_lo <= 0 & a.ci_hi >= 0), ...
        sum(f.ci_hi < 0), sum(f.ci_lo > 0), sum(f.ci_lo <= 0 & f.ci_hi >= 0), sum(f.n_samples == 0), ...
        sum(sign(a.slope(ia(both))) == sign(f.slope(ib(both)))), sum(both));
end
fprintf('\nEnd %s\n', datestr(now)); diary off

function s = fl(q)
    s = ' ';
    if isempty(q) || ~isfinite(q.slope), s = '?'; return; end
    if q.ci_hi < 0, s = '-'; elseif q.ci_lo > 0, s = '+'; end
end
