% Used for: interannual GLMM by biome on the Fig. S8a samples (Fig. R4)
% Original: Manuscript\MLD\manuscript\nature\NCC\R2\code for revision\interannual_biome\interannual_biome.m (data paths updated to I:\Prj; computation unchanged)
% interannual_biome.m
% Biome version of the interannual analysis (HZ 2026-11-05, see revision_log_R2.md), on the Fig. S8a samples:
% grid cell x calendar month pairs with BOTH NPP and MLD-NEMO Mann-Kendall p < 0.05, 14 non-polar biomes,
% hemisphere-aware seasons, exactly as R1\code for revision\mld_npp_relationship_biome\compute_biome_sensitivity.m.
% Per sample: 20 monthly values (2003-2022) of ln NPP and ln MLD, linear detrend, OLS slope.
% GLMM per NPP model x season: lnNPP ~ lnMLD + (1 + lnMLD | biome), biome slope = fixed + biome random effect.
% CbPM, VGPM, CAFE. READ ONLY on inputs. Haipeng Zhao / Claude, 2026-11-05 (R2)

clear; clc;
OUT = fileparts(mfilename('fullpath'));
TR = 'pct';                                   % 'log' (first run) or 'pct' (percent anomaly, HZ 2026-11-05)
sfx = ''; if strcmp(TR, 'pct'), sfx = '_pct'; end
tf = @(v) log(v); if strcmp(TR, 'pct'), tf = @(v) (v - mean(v)) / mean(v) * 100; end
diary(fullfile(OUT, ['interannual_biome' sfx '_log.txt'])); diary on
fprintf('Start %s, transform %s\n', datestr(now), TR);
R1B = fullfile(fileparts(fileparts(OUT)), 'Extended_data_figures', 'FigS8_biome_sensitivity');
pthr = 0.05; YEARS = 2003:2022; NY = 20;
mmN = [3 4 5; 6 7 8; 9 10 11; 12 1 2];  mmS = [9 10 11; 12 1 2; 3 4 5; 6 7 8];
seasName = {'Spring','Summer','Autumn','Winter'};
biomes = {'A EQU','EP EQU','WP EQU','NA STPS','SA STPS','NP STPS','SP STPS','IND STPS', ...
          'NA STSS','NP STSS','SO STSS','NA SPSS','NP SPSS','SO SPSS'};       % same list and order as R1
NB = numel(biomes);

Tr = load(fullfile(R1B, 'grid_monthly_trends_mld_npp_nemo.mat'));             % CbPM + MLD trends (R1)
Tv = load(fullfile(OUT, 'grid_monthly_trends_npp_vgpm_cafe.mat'));             % VGPM + CAFE trends (R2)
S8 = load(fullfile(R1B, 'senMLD_biome_nemo_sig005.mat'), 'senMLD', 'biomes');
assert(isequal(S8.biomes, biomes), 'biome list differs from S8a');
lat = double(Tr.lat(:)); isN = lat >= 0;
trend.CbPM = {Tr.npp_trend_pct, Tr.npp_trend_p};
trend.VGPM = {Tv.vgpm_trend_pct, Tv.vgpm_trend_p};
trend.CAFE = {Tv.cafe_trend_pct, Tv.cafe_trend_p};
src = {'CbPM', 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\grid_npp_nemo\npp_YYYY.mat';
       'VGPM', 'I:\Prj\npp_vgpm_1M_1dg\grid_npp_vgpm_1M_1dg_mat\npp_vgpm_YYYY.mat';
       'CAFE', 'I:\Prj\npp_cafe_1M_1dg\grid_npp_nemo_1M_1dg_mat\npp_YYYY.mat'};

% biome cells as (lat,lon) linear indices (idx is into a (360,180) = [lon,lat] array)
bcell = cell(NB,1);
for b = 1:NB
    d = load(fullfile('I:\Prj\MLD_1M_NEMO\biomes_1M_1dg\mld\2003', ['mld_' biomes{b} '.mat']), 'idx');
    [ilon, ilat] = ind2sub([360 180], round(d.idx(:))); bcell{b} = sub2ind([180 360], ilat, ilon);
end

MLDm = nan(180,360,12,NY,'single');
for y = 1:NY, d = load(['I:\Prj\MLD_1M_NEMO\one_d03\mld\mld_' num2str(YEARS(y)) '.mat'], 'mldG1deg'); MLDm(:,:,:,y) = d.mldG1deg; end
MLDm = reshape(MLDm, 180*360, 12, NY);

R = {}; SAMP = struct();
for k = 1:3
    NPPm = nan(180,360,12,NY,'single');
    for y = 1:NY, d = load(strrep(src{k,2}, 'YYYY', num2str(YEARS(y))), 'nppG1dg'); NPPm(:,:,:,y) = d.nppG1dg; end
    NPPm = reshape(NPPm, 180*360, 12, NY);
    tp = trend.(src{k,1}){1}; pp = trend.(src{k,1}){2};
    for s = 1:4
        rows = []; info = [];                       % info: biome, cell, month, OLS slope
        for b = 1:NB
            c = bcell{b}; [ilat, ~] = ind2sub([180 360], c); nC = 0; nBad = 0;
            for mi = 1:3
                for h = 1:2
                    if h == 1, sel = isN(ilat); m = mmN(s,mi); else, sel = ~isN(ilat); m = mmS(s,mi); end
                    cc = c(sel);
                    nt = double(tp(:,:,m)); np = double(pp(:,:,m)); mt = double(Tr.mld_trend_pct(:,:,m)); mp = double(Tr.mld_trend_p(:,:,m));
                    ok = ~isnan(nt(cc)) & np(cc) < pthr & mp(cc) < pthr;
                    ratio = nt(cc) ./ mt(cc); ok = ok & isfinite(ratio);          % as in R1: drop inf / NaN ratios
                    cc = cc(ok); nC = nC + numel(cc);
                    for q = 1:numel(cc)
                        yv = double(squeeze(NPPm(cc(q), m, :))); xv = double(squeeze(MLDm(cc(q), m, :)));
                        if any(~isfinite([yv; xv])) || any(xv <= 0) || any(yv < 0) || (strcmp(TR, 'log') && any(yv <= 0)), nBad = nBad + 1; continue; end
                        yd = detrend(tf(yv)); xd = detrend(tf(xv));
                        xc = xd - mean(xd); bq = sum(xc .* (yd - mean(yd))) / sum(xc.^2);
                        info(end+1,:) = [b cc(q) m bq]; %#ok<AGROW>
                        rows = [rows; yd xd repmat([b size(info,1)], NY, 1)]; %#ok<AGROW>
                    end
                end
            end
            if k == 1 && nC ~= S8.senMLD(b,s,3) && ~(nC == 0 && isnan(S8.senMLD(b,s,3)))
                fprintf('  COUNT MISMATCH CbPM %s %s: %d vs S8a %d\n', biomes{b}, seasName{s}, nC, S8.senMLD(b,s,3));
            end
            if nBad > 0, fprintf('  %s %s %s: %d samples dropped (non-positive or missing values)\n', src{k,1}, biomes{b}, seasName{s}, nBad); end
        end
        SAMP.(src{k,1}){s} = info;
        T = table(rows(:,1), rows(:,2), categorical(rows(:,3)), 'VariableNames', {'lnNPP','lnMLD','biome'});
        out = nan(NB, 3); fx = nan(1,3);
        try
            mdl = fitlme(T, 'lnNPP ~ lnMLD + (1 + lnMLD | biome)'); ci = coefCI(mdl);
            b0 = mdl.Coefficients.Estimate(2); se0 = mdl.Coefficients.SE(2); fx = [b0 ci(2,:)];
            [~, ~, st] = randomEffects(mdl); isB = strcmp(st.Name, 'lnMLD');
            lev = str2double(st.Level(isB)); bb = st.Estimate(isB); sp = st.SEPred(isB);
            for j = 1:numel(lev), e = b0 + bb(j); hw = 1.96 * sqrt(se0^2 + sp(j)^2); out(lev(j),:) = [e e - hw e + hw]; end
        catch e, fprintf('  LME fail %s %s: %s\n', src{k,1}, seasName{s}, e.message);
        end
        R(end+1,:) = {src{k,1}, seasName{s}, 'All biomes', size(info,1), mean(info(:,4) < 0), median(info(:,4)), fx(1), fx(2), fx(3)}; %#ok<SAGROW>
        for b = 1:NB
            v = info(info(:,1) == b, 4);
            R(end+1,:) = {src{k,1}, seasName{s}, biomes{b}, numel(v), mean(v < 0), median(v), out(b,1), out(b,2), out(b,3)}; %#ok<SAGROW>
        end
        fprintf('%s %s: %d samples\n', src{k,1}, seasName{s}, size(info,1));
    end
end
Tab = cell2table(R, 'VariableNames', {'npp','season','biome','n','pneg_ols','median_ols','slope','ci_lo','ci_hi'});
Tab.slope(Tab.n == 0) = NaN; Tab.ci_lo(Tab.n == 0) = NaN; Tab.ci_hi(Tab.n == 0) = NaN;
writetable(Tab, fullfile(OUT, ['interannual_biome' sfx '.csv']));
save(fullfile(OUT, ['interannual_biome' sfx '.mat']), 'Tab', 'SAMP', 'biomes', 'seasName');

ordB = {'NA SPSS','NP SPSS','NA STSS','NP STSS','NA STPS','NP STPS','EP EQU','A EQU','WP EQU','IND STPS','SA STPS','SP STPS','SO STSS','SO SPSS'};
fprintf('\nGLMM biome slope [95%% CI] (n samples) and %% samples negative (OLS); biomes north to south, cols Sp Su Au Wi\n');
for k = 1:3
    fprintf('-- %s\n', src{k,1});
    for un = [{'All biomes'}, ordB]
        fprintf('  %-10s', un{1});
        for s = 1:4
            q = Tab(strcmp(Tab.npp, src{k,1}) & strcmp(Tab.biome, un{1}) & strcmp(Tab.season, seasName{s}), :);
            fl = ' '; if q.ci_hi < 0, fl = '-'; elseif q.ci_lo > 0, fl = '+'; end
            fprintf('  %5.2f[%5.2f,%5.2f]%s(%4d,%3.0f%%)', q.slope, q.ci_lo, q.ci_hi, fl, q.n, 100*q.pneg_ols);
        end
        fprintf('\n');
    end
end
fprintf('\nBiome-seasons (56): CI<0 / CI>0 / spans 0 / no samples; pooled samples negative\n');
for k = 1:3
    t = Tab(strcmp(Tab.npp, src{k,1}) & ~strcmp(Tab.biome, 'All biomes'), :);
    fprintf('  %s: %2d / %2d / %2d / %2d ; %.0f%% of %d samples negative\n', src{k,1}, sum(t.ci_hi < 0), sum(t.ci_lo > 0), ...
        sum(t.ci_lo <= 0 & t.ci_hi >= 0), sum(t.n == 0), 100*sum(t.pneg_ols .* t.n, 'omitnan')/sum(t.n), sum(t.n));
end
fprintf('\nEnd %s\n', datestr(now)); diary off
