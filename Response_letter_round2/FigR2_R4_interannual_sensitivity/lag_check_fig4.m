% Used for: lag check, lag_check_fig4.mat (read by glmm_singular_all_figures.m; not shown in the letter)
% Original: Manuscript\MLD\manuscript\nature\NCC\R2\code for revision\interannual_fig4\lag_check_fig4.m (data paths updated to I:\Prj; computation unchanged)
% lag_check_fig4.m
% Lag effect in the interannual NPP-MLD relationship (pre-specified in revision_log_R2.md, 2026-11-06).
% Same Fig. 4 samples (lat row x month m, paired significant trends p < 0.1), zonal monthly series, CbPM / VGPM / CAFE.
% NPP in month m of year y vs MLD at: lag 0 (month m), 1, 2, 3 months earlier (continuous monthly axis),
% the mean of months m-3..m-1 (one-season lag), and month m of year y-1 (one-year lag).
% Each series: percent anomaly over the years used, linear detrend, per-sample OLS slope.
% GLMM per model x season x lag: y ~ x + (1 + x | band), band slope = fixed + band effect, approx 95% CI.
% Per sample also y ~ x(lag 0) + x(season lag): lag coefficient and sum of both.
% READ ONLY on inputs. Haipeng Zhao / Claude, 2026-11-06 (R2)

clear; clc;
OUT = fileparts(mfilename('fullpath'));
diary(fullfile(OUT, 'lag_check_fig4_log.txt')); diary on
fprintf('Start %s\n', datestr(now));
load(fullfile(OUT, 'interannual_fig4_samples_pct.mat'), 'S');
mmN = [3 4 5; 6 7 8; 9 10 11; 12 1 2];  mmS = [9 10 11; 12 1 2; 3 4 5; 6 7 8];
seasName = {'Spring','Summer','Autumn','Winter'};
latBand = [60 40; 40 20; 20 0; 0 -20; -20 -40; -40 -60];
bandName = {'40-60N','20-40N','0-20N','0-20S','20-40S','40-60S'};
models = {'CbPM','VGPM','CAFE'};
nppFile = {'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\nppLat_monthly_nemo.mat', 'I:\Prj\npp_vgpm_1M_1dg\nppLat_monthly_vgpm.mat', ...
           'I:\Prj\npp_cafe_1M_1dg\nppLat_monthly_nemo.mat'};
lagName = {'lag0','lag1m','lag2m','lag3m','season','lag1y'};
lagIdx = {0, 1, 2, 3, 1:3, 12};                       % months back on the continuous axis (mean if several)
NL = numel(lagName); NY = 20;
pa = @(v) (v - mean(v)) / mean(v) * 100;

dm = load('I:\Prj\MLD_1M_NEMO\mldLat_monthly_nemo.mat');
latS = double(dm.lat(:)); keepM = latS < 60 & latS > -60; mldLat = double(dm.mldLat(keepM, :)); mmM = dm.mm(:)';
assert(isequal(mmM, repmat(1:12, 1, NY)), 'MLD month axis is not continuous Jan 2003 - Dec 2022');
R = {}; RS = {}; SL = struct();
for k = 1:3
    dn = load(nppFile{k}); latN = double(dn.lat(:)); keep = latN < 60 & latN > -60;
    nppLat = double(dn.nppLat(keep, 1:NY*12)); nLat = sum(keep);
    assert(isequal(dn.mm(1:NY*12), mmM) || isequal(dn.mm(1:NY*12)', mmM), 'NPP month axis mismatch');
    ok = S.(models{k}).ok; latF = S.(models{k}).latF;
    for s = 1:4
        rows = cell(NL, 1); bS = cell(NL, 6); both = cell(6, 1);
        for u = 1:6
            [r, c] = sampleIdx(ok, latF, latBand(u,:), mmN(s,:), mmS(s,:));
            for q = 1:numel(r)
                rS = nLat + 1 - r(q); i0 = ((1:NY) - 1) * 12 + c(q);
                X = nan(NY, NL);
                for L = 1:NL
                    for y = 1:NY
                        ii = i0(y) - lagIdx{L};
                        if all(ii >= 1), X(y, L) = mean(mldLat(rS, ii)); end
                    end
                end
                yv = nppLat(rS, i0)';
                if any(~isfinite(yv)) || any(yv <= 0) || any(X(:,1) <= 0), continue; end
                for L = 1:NL
                    v = isfinite(X(:, L));
                    yd = detrend(pa(yv(v))); xd = detrend(pa(X(v, L)));
                    bS{L, u}(end+1) = olsb(yd, xd);
                    rows{L} = [rows{L}; yd xd repmat(u, sum(v), 1)];
                end
                v = isfinite(X(:, 5));                                  % lag 0 + season lag, same years
                yd = detrend(pa(yv(v))); x0 = detrend(pa(X(v, 1))); xs = detrend(pa(X(v, 5)));
                cf = [x0 xs ones(sum(v), 1)] \ yd; both{u}(end+1, :) = [cf(2) cf(1) + cf(2)];
            end
        end
        SL.(models{k}){s} = bS;                            % per-sample slopes, lag x band
        for L = 1:NL
            out = nan(6, 3);
            if ~isempty(rows{L})
                T = table(rows{L}(:,1), rows{L}(:,2), categorical(rows{L}(:,3)), 'VariableNames', {'y','x','band'});
                try
                    m = fitlme(T, 'y ~ x + (1 + x | band)'); b0 = m.Coefficients.Estimate(2); se0 = m.Coefficients.SE(2);
                    [~, ~, st] = randomEffects(m); isB = strcmp(st.Name, 'x');
                    lev = str2double(st.Level(isB)); bb = st.Estimate(isB); sp = st.SEPred(isB);
                    for j = 1:numel(lev), e = b0 + bb(j); hw = 1.96 * sqrt(se0^2 + sp(j)^2); out(lev(j), :) = [e e - hw e + hw]; end
                catch e, fprintf('  LME fail %s %s %s: %s\n', models{k}, seasName{s}, lagName{L}, e.message);
                end
            end
            for u = 1:6
                if isempty(bS{L, u}), continue; end
                R(end+1,:) = {models{k}, bandName{u}, seasName{s}, lagName{L}, numel(bS{L, u}), median(bS{L, u}), out(u,1), out(u,2), out(u,3)}; %#ok<SAGROW>
            end
        end
        for u = 1:6
            if isempty(both{u}), continue; end
            RS(end+1,:) = {models{k}, bandName{u}, seasName{s}, size(both{u}, 1), median(both{u}(:,1)), mean(both{u}(:,1) < 0), ...
                median(both{u}(:,2)), mean(both{u}(:,2) < 0)}; %#ok<SAGROW>
        end
    end
end
Tab = cell2table(R, 'VariableNames', {'npp','band','season','lag','n','ols_median','glmm','glmm_lo','glmm_hi'});
TabB = cell2table(RS, 'VariableNames', {'npp','band','season','n','lagcoef_median','lagcoef_pneg','sum_median','sum_pneg'});
writetable(Tab, fullfile(OUT, 'lag_check_fig4.csv')); writetable(TabB, fullfile(OUT, 'lag_check_fig4_both.csv'));
save(fullfile(OUT, 'lag_check_fig4.mat'), 'Tab', 'TabB', 'SL', 'lagName');

fprintf('\nGLMM band slope by lag (flag -/+ if CI excludes 0) | lag0+season model: lag coef, sum (medians)\n');
fprintf('  %-7s %-7s %4s', 'band', 'season', 'n'); fprintf(' %9s', lagName{:}); fprintf(' | %7s %7s\n', 'lagcoef', 'sum');
for k = 1:3
    fprintf('-- %s\n', models{k});
    for u = 1:6
        for s = 1:4
            q = Tab(strcmp(Tab.npp, models{k}) & strcmp(Tab.band, bandName{u}) & strcmp(Tab.season, seasName{s}), :);
            if isempty(q), continue; end
            fprintf('  %-7s %-7s %4d', bandName{u}, seasName{s}, q.n(1));
            for L = 1:NL
                z = q(strcmp(q.lag, lagName{L}), :); fl = ' ';
                if z.glmm_hi < 0, fl = '-'; elseif z.glmm_lo > 0, fl = '+'; end
                fprintf(' %8.2f%s', z.glmm, fl);
            end
            b = TabB(strcmp(TabB.npp, models{k}) & strcmp(TabB.band, bandName{u}) & strcmp(TabB.season, seasName{s}), :);
            fprintf(' | %7.2f %7.2f\n', b.lagcoef_median, b.sum_median);
        end
    end
end
fprintf('\nBand-seasons GLMM CI<0 / CI>0 / spans 0, per lag\n');
for k = 1:3
    fprintf('  %s:', models{k});
    for L = 1:NL
        t = Tab(strcmp(Tab.npp, models{k}) & strcmp(Tab.lag, lagName{L}), :);
        fprintf('  %s %d/%d/%d', lagName{L}, sum(t.glmm_hi < 0), sum(t.glmm_lo > 0), sum(t.glmm_lo <= 0 & t.glmm_hi >= 0));
    end
    fprintf('\n');
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
function b = olsb(y, x)
    xc = x - mean(x); b = sum(xc .* (y - mean(y))) / sum(xc.^2);
end
