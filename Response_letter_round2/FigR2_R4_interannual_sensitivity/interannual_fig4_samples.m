% Used for: interannual (detrended percent anomaly) samples on the Fig. 4 samples (Fig. R2)
% Original: Manuscript\MLD\manuscript\nature\NCC\R2\code for revision\interannual_fig4\interannual_fig4_samples.m (data paths updated to I:\Prj; computation unchanged)
% interannual_fig4_samples.m
% Reviewer 2 comment 1 (HZ decision 2026-11-05, see revision_log_R2.md): interannual NPP-MLD relationship on
% EXACTLY the samples behind Fig. 4: latitude row x calendar month (zonal means, 60S-60N) with concurrent
% significant NPP and MLD-NEMO trends (Mann-Kendall p < 0.1). CbPM, VGPM, CAFE separately.
% Per sample: ln NPP and ln MLD, 20 yearly values, linear detrend, OLS slope (elasticity).
% Per band x season: n, % negative, median slope; mixed model lnNPP ~ lnMLD + (1 + lnMLD | sample) where n >= 3.
% Absolute-value slopes kept for reference only. READ ONLY on inputs. Haipeng Zhao / Claude, 2026-11-05 (R2)

clear; clc;
OUT = fileparts(mfilename('fullpath'));
ROOT = 'D:\Users\zhp19\OneDrive - Indiana State University\Research';
if ~isfolder(ROOT), ROOT = 'E:\OneDrive - Indiana State University\Research'; end      % laptop
PRJ = fullfile(ROOT, 'Projects', '3D bio model');
addpath(fullfile(ROOT, 'MATLAB toolbox', 'ktaub'));
pthr = 0.1;
TR = 'pct';                                   % 'log' (first run) or 'pct' (percent anomaly, HZ 2026-11-05)
sfx = ''; if strcmp(TR, 'pct'), sfx = '_pct'; end
tf = @(v) log(v); if strcmp(TR, 'pct'), tf = @(v) (v - mean(v)) / mean(v) * 100; end
diary(fullfile(OUT, ['interannual_fig4_samples' sfx '_log.txt'])); diary on
fprintf('Start %s, transform %s\n', datestr(now), TR);
mmN = [3 4 5; 6 7 8; 9 10 11; 12 1 2];  mmS = [9 10 11; 12 1 2; 3 4 5; 6 7 8];
seasName = {'Spring','Summer','Autumn','Winter'};
latBand = [60 40; 40 20; 20 0; 0 -20; -20 -40; -40 -60];
bandName = {'40-60N','20-40N','0-20N','0-20S','20-40S','40-60S'};
src = {'CbPM', fullfile(fileparts(fileparts(OUT)), '4_trends_and_sensitivity', 'NPPtrends.mat'),                 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\nppLat_monthly_nemo.mat';
       'VGPM', 'I:\Prj\npp_vgpm_1M_1dg\VGPM_NPPtrends.mat', ...
                                                                       'I:\Prj\npp_vgpm_1M_1dg\nppLat_monthly_vgpm.mat';
       'CAFE', fullfile(fileparts(fileparts(OUT)), 'Supplementary_information', 'SuppFig6_CAFE_month_latitude_trends', 'CAFE_NPPtrends.mat'), 'I:\Prj\npp_cafe_1M_1dg\nppLat_monthly_nemo.mat'};
mldFile = 'I:\Prj\MLD_1M_NEMO\mldLat_monthly_nemo.mat';
haveI = isfile(mldFile);
if ~haveI, fprintf('I:\\Prj not available: sample counts only.\n'); end

%% 1. Fig. 4 samples (lat row x month with paired significant trends), as in Code\senMLD2_ot.m
S = struct();
for k = 1:3
    T = load(src{k,2}, 'lat', 'npp_trends_sig_nemo', 'mld_trends_sig_nemo', 'npp_trends_pct_nemo', 'mld_trends_pct_nemo');
    latF = flip(double(T.lat(:)));                         % trend files: row 1 = north (rot90 in the trend script)
    ok = T.npp_trends_sig_nemo < pthr & T.mld_trends_sig_nemo < pthr & isfinite(T.npp_trends_pct_nemo);
    S.(src{k,1}).ok = ok; S.(src{k,1}).latF = latF; S.(src{k,1}).T = T;
end
fprintf('\nFig. 4 sample counts (lat row x month, paired significant trends p<%.2f), rows 40-60N ... 40-60S, cols Sp Su Au Wi\n', pthr);
for k = 1:3
    fprintf('-- %s\n', src{k,1}); ok = S.(src{k,1}).ok; latF = S.(src{k,1}).latF;
    for u = 1:6
        fprintf('  %-7s', bandName{u});
        for s = 1:4
            [r, ~] = sampleIdx(ok, latF, latBand(u,:), mmN(s,:), mmS(s,:));
            fprintf(' %3d', numel(r));
        end
        fprintf('\n');
    end
end
if ~haveI, fprintf('\nEnd %s\n', datestr(now)); diary off; return; end

%% 2. Interannual (detrended) regression on those samples
dm = load(mldFile);
latS = double(dm.lat(:)); keepM = latS < 60 & latS > -60; mldLat = double(dm.mldLat(keepM, :)); latS = latS(keepM);
mmM = dm.mm(:)';
R = {}; P = struct();
for k = 1:3
    dn = load(src{k,3});
    latN = double(dn.lat(:)); keep = latN < 60 & latN > -60; nppLat = double(dn.nppLat(keep, :)); latN = latN(keep);
    if numel(dn.mm) > numel(mmM)            % VGPM zonal file is 2003-2023: keep 2003-2022, as in R1 npp_mld_trends_vgpm_ot.m
        nppLat = nppLat(:, 1:numel(mmM)); dn.mm = dn.mm(1:numel(mmM));
    end
    assert(isequal(latN, latS), 'lat mismatch between NPP and MLD series');
    assert(isequal(dn.mm(:)', mmM), 'month axis mismatch');
    latF = S.(src{k,1}).latF; assert(max(abs(flip(latN) - latF)) < 1e-3, 'trend-file lat is not the flipped series lat');
    ok = S.(src{k,1}).ok; T = S.(src{k,1}).T;
    nLat = numel(latN); bD = nan(nLat, 12); bA = bD; pD = bD; dP = nan(nLat, 12, 2);
    for rF = 1:nLat
        rS = nLat + 1 - rF;                                % series row for trend-file row rF
        for j = 1:12
            if ~ok(rF, j), continue; end
            y = nppLat(rS, mmM == j)'; x = mldLat(rS, mmM == j)';
            if any(~isfinite([y; x])) || any([y; x] <= 0), continue; end
            t = (1:numel(y))';
            [~,~,~,pn] = ktaub([t y], 0.9999, 0); [~,~,~,pm] = ktaub([t x], 0.9999, 0);
            dP(rF, j, :) = [pn - T.npp_trends_sig_nemo(rF, j), pm - T.mld_trends_sig_nemo(rF, j)];
            ly = tf(y); lx = tf(x);
            bA(rF, j) = olsb(ly, lx);
            yd = detrend(ly); xd = detrend(lx);
            [bD(rF, j), pD(rF, j)] = olsbp(yd, xd);
            P.(src{k,1}).series{rF, j} = [yd xd];
        end
    end
    fprintf('\n%s: ground-truth check, max |recomputed MK p - stored p|: NPP %.2g, MLD %.2g\n', src{k,1}, max(abs(dP(:,:,1)), [], 'all'), max(abs(dP(:,:,2)), [], 'all'));
    P.(src{k,1}).bD = bD; P.(src{k,1}).bA = bA; P.(src{k,1}).pD = pD;
    for u = 0:6
        for s = 1:4
            if u == 0, [r, c] = sampleIdx(ok, latF, [60 -60], mmN(s,:), mmS(s,:)); un = 'All 60S-60N';
            else, [r, c] = sampleIdx(ok, latF, latBand(u,:), mmN(s,:), mmS(s,:)); un = bandName{u}; end
            li = sub2ind(size(bD), r, c); bd = bD(li); ba = bA(li); pd = pD(li); v = isfinite(bd);
            out = nan(1,3);
            if sum(v) >= 3
                rows = []; vi = find(v);
                for q = 1:numel(vi), z = P.(src{k,1}).series{r(vi(q)), c(vi(q))}; rows = [rows; z repmat(q, size(z,1), 1)]; end %#ok<AGROW>
                Tb = table(rows(:,1), rows(:,2), categorical(rows(:,3)), 'VariableNames', {'lnNPP','lnMLD','sample'});
                try
                    m = fitlme(Tb, 'lnNPP ~ lnMLD + (1 + lnMLD | sample)'); ci = coefCI(m);
                    out = [m.Coefficients.Estimate(2) ci(2,:)];
                catch e, fprintf('  LME fail %s %s %s: %s\n', src{k,1}, un, seasName{s}, e.message); end
            end
            R(end+1,:) = [{src{k,1}, un, seasName{s}, sum(v), mean(bd(v) < 0), mean(bd(v) < 0 & pd(v) < pthr), mean(bd(v) > 0 & pd(v) < pthr), ...
                median(bd(v)), mean(ba(v) < 0), median(ba(v))}, num2cell(out)]; %#ok<SAGROW>
        end
    end
end
Tab = cell2table(R, 'VariableNames', {'npp','region','season','n','pneg_detr','psig_neg_detr','psig_pos_detr','median_detr', ...
    'pneg_abs','median_abs','lme_detr','lme_detr_lo','lme_detr_hi'});
writetable(Tab, fullfile(OUT, ['interannual_fig4_samples' sfx '.csv']));
save(fullfile(OUT, ['interannual_fig4_samples' sfx '.mat']), 'Tab', 'P', 'S');

fprintf('\nDETRENDED: %% samples negative [median slope] (n), rows = region, cols Sp Su Au Wi\n');
for k = 1:3
    fprintf('-- %s\n', src{k,1});
    for un = [{'All 60S-60N'}, bandName]
        fprintf('  %-12s', un{1});
        for s = 1:4
            q = Tab(strcmp(Tab.npp, src{k,1}) & strcmp(Tab.region, un{1}) & strcmp(Tab.season, seasName{s}), :);
            fprintf(' %4.0f%% [%5.2f] (%3d)', 100*q.pneg_detr, q.median_detr, q.n);
        end
        fprintf('\n');
    end
end
fprintf('\nBand x season (24): median detrended slope <0 / >0 / no samples; LME CI<0 / CI>0 / spans 0 / not fitted\n');
for k = 1:3
    t = Tab(strcmp(Tab.npp, src{k,1}) & ~strcmp(Tab.region, 'All 60S-60N'), :);
    fprintf('  %s: %2d / %2d / %2d ; LME %2d / %2d / %2d / %2d ; pooled samples negative %.0f%% (n=%d)\n', src{k,1}, ...
        sum(t.median_detr < 0), sum(t.median_detr > 0), sum(t.n == 0), sum(t.lme_detr_hi < 0), sum(t.lme_detr_lo > 0), ...
        sum(t.lme_detr_lo <= 0 & t.lme_detr_hi >= 0), sum(isnan(t.lme_detr)), 100*sum(t.pneg_detr .* t.n, 'omitnan')/sum(t.n), sum(t.n));
end
fprintf('\nEnd %s\n', datestr(now)); diary off

function [r, c] = sampleIdx(ok, latF, band, mN, mS)
    % rows/cols of paired-significant samples in a band and season (hemisphere-aware months), as in senMLD2_ot.m
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
function [b, p] = olsbp(y, x)
    n = numel(y); xc = x - mean(x); yc = y - mean(y); sxx = sum(xc.^2);
    b = sum(xc .* yc) / sxx; res = yc - b * xc; se = sqrt(sum(res.^2) / (n - 2) / sxx);
    p = 2 * tcdf(-abs(b / se), n - 2);
end
