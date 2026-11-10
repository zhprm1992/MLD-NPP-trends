% Used for: trend-based sensitivity of the same samples for CbPM, VGPM, CAFE (Fig. R2 top row)
% Original: Manuscript\MLD\manuscript\nature\NCC\R2\code for revision\interannual_fig4\trend_vs_interannual.m (data paths updated to I:\Prj; computation unchanged)
% trend_vs_interannual.m
% Same Fig. 4 samples (lat row x month, paired significant trends p < 0.1): trend-to-trend sensitivity vs
% interannual (detrended percent-anomaly) slope, per NPP model. Pre-specified in revision_log_R2.md (2026-11-06).
% Haipeng Zhao / Claude, 2026-11-06 (R2)

clear; clc;
OUT = fileparts(mfilename('fullpath'));
diary(fullfile(OUT, 'trend_vs_interannual_log.txt')); diary on
fprintf('Start %s\n', datestr(now));
load(fullfile(OUT, 'interannual_fig4_samples_pct.mat'), 'P', 'S');
G = load(fullfile(OUT, 'glmm_band_group_pct.mat'), 'Tab'); G = G.Tab(strcmp(G.Tab.model, 'main'), :);
mmN = [3 4 5; 6 7 8; 9 10 11; 12 1 2];  mmS = [9 10 11; 12 1 2; 3 4 5; 6 7 8];
seasName = {'Spring','Summer','Autumn','Winter'};
latBand = [60 40; 40 20; 20 0; 0 -20; -20 -40; -40 -60];
bandName = {'40-60N','20-40N','0-20N','0-20S','20-40S','40-60S'};
models = {'CbPM','VGPM','CAFE'};
R = {}; ALL = struct();
for k = 1:3
    T = S.(models{k}).T; ok = S.(models{k}).ok; latF = S.(models{k}).latF;
    A = zeros(0, 6);                               % per sample: band, season, trendSens, iaSlope, ratioMLD, ratioNPP
    for u = 1:6
        for s = 1:4
            [r, c] = sampleIdx(ok, latF, latBand(u,:), mmN(s,:), mmS(s,:));
            for q = 1:numel(r)
                z = P.(models{k}).series{r(q), c(q)};
                if isempty(z), continue; end
                tn = T.npp_trends_pct_nemo(r(q), c(q)); tm = T.mld_trends_pct_nemo(r(q), c(q));
                rM = std(z(:,2)) / (abs(tm) * 19); rN = std(z(:,1)) / (abs(tn) * 19);
                A(end+1,:) = [u s tn/tm P.(models{k}).bD(r(q), c(q)) rM rN]; %#ok<AGROW>
            end
            sel = A(:,1) == u & A(:,2) == s; g = G(strcmp(G.npp, models{k}) & strcmp(G.band, bandName{u}) & strcmp(G.season, seasName{s}), :);
            if any(sel)
                R(end+1,:) = {models{k}, bandName{u}, seasName{s}, sum(sel), median(A(sel,3)), g.slope, median(A(sel,4)), ...
                    mean(A(sel,3) < 0 & A(sel,4) < 0), mean(A(sel,3) < 0 & A(sel,4) > 0), mean(A(sel,3) > 0 & A(sel,4) > 0), mean(A(sel,3) > 0 & A(sel,4) < 0)}; %#ok<AGROW>
            end
        end
    end
    ALL.(models{k}) = A;
end
Tab = cell2table(R, 'VariableNames', {'npp','band','season','n','median_trend_sens','glmm_ia_slope','median_ia_slope', ...
    'f_trendNeg_iaNeg','f_trendNeg_iaPos','f_trendPos_iaPos','f_trendPos_iaNeg'});
writetable(Tab, fullfile(OUT, 'trend_vs_interannual.csv'));
save(fullfile(OUT, 'trend_vs_interannual.mat'), 'Tab', 'ALL');

fprintf('\nPer band x season: median trend sensitivity | GLMM interannual slope (n)\n');
for k = 1:3
    fprintf('-- %s\n', models{k});
    for u = 1:6
        fprintf('  %-7s', bandName{u});
        for s = 1:4
            q = Tab(strcmp(Tab.npp, models{k}) & strcmp(Tab.band, bandName{u}) & strcmp(Tab.season, seasName{s}), :);
            if isempty(q), fprintf('  %-22s', '   -'); else, fprintf('  %6.2f | %5.2f (%2d)  ', q.median_trend_sens, q.glmm_ia_slope, q.n); end
        end
        fprintf('\n');
    end
end
fprintf('\nPer sample sign table (trend / interannual): -/-  -/+  +/+  +/-   (n)\n');
for k = 1:3
    A = ALL.(models{k}); a = A(:,3); b = A(:,4);
    fprintf('  %s: %3d %3d %3d %3d  (%d)\n', models{k}, sum(a<0&b<0), sum(a<0&b>0), sum(a>0&b>0), sum(a>0&b<0), numel(a));
end
fprintf('\nVariability test: interannual SD / total trend change (median [IQR]), MLD and NPP\n');
for k = 1:3
    A = ALL.(models{k}); a = A(:,3); b = A(:,4);
    g1 = a < 0 & b < 0; g2 = a < 0 & b > 0;
    for v = 5:6
        if v == 5, vn = 'MLD'; else, vn = 'NPP'; end
        p = NaN; if sum(g1) > 2 && sum(g2) > 2, p = ranksum(A(g1,v), A(g2,v)); end
        fprintf('  %s %s: trend-/ia- %.2f [%.2f %.2f] (n=%d) | trend-/ia+ %.2f [%.2f %.2f] (n=%d) | rank-sum p=%.3f\n', models{k}, vn, ...
            median(A(g1,v)), prctile(A(g1,v),25), prctile(A(g1,v),75), sum(g1), median(A(g2,v)), prctile(A(g2,v),25), prctile(A(g2,v),75), sum(g2), p);
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
