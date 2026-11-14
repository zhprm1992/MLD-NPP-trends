% Used for: biome interannual GLMM on all grid cells and months, no significance filter (Fig. R4 bottom row)
% Original: Manuscript\MLD\manuscript\nature\NCC\R2\code for revision\interannual_biome\interannual_biome_allcells.m (data paths updated to I:\Prj; computation unchanged)
% interannual_biome_allcells.m (2026-11-10), internal check for the meeting (see revision_log_R2.md).
% Same inputs, biome cells, seasons, transform and model as interannual_biome.m (Figure R4, TR = 'pct'),
% with ONE change: no significance filter. Every biome grid cell x calendar month with 20 finite values,
% MLD > 0, NPP >= 0 and mean NPP > 0 is a sample. Vectorized sample building (many more samples).
% Model per NPP model x season: lnNPP ~ lnMLD + (1 + lnMLD | biome) on detrended percent anomalies
% (variable names kept from interannual_biome.m). READ ONLY on inputs. Haipeng Zhao / Claude, 2026-11-10.

clear; clc;
OUT = fileparts(mfilename('fullpath'));
diary(fullfile(OUT, 'interannual_biome_allcells_log.txt')); diary on
fprintf('Start %s\n', datestr(now));
YEARS = 2003:2022; NY = 20;
mmN = [3 4 5; 6 7 8; 9 10 11; 12 1 2];  mmS = [9 10 11; 12 1 2; 3 4 5; 6 7 8];
seasName = {'Spring','Summer','Autumn','Winter'};
biomes = {'A EQU','EP EQU','WP EQU','NA STPS','SA STPS','NP STPS','SP STPS','IND STPS', ...
          'NA STSS','NP STSS','SO STSS','NA SPSS','NP SPSS','SO SPSS'};
NB = numel(biomes);
src = {'CbPM', 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\grid_npp_nemo\npp_YYYY.mat';
       'VGPM', 'I:\Prj\npp_vgpm_1M_1dg\grid_npp_vgpm_1M_1dg_mat\npp_vgpm_YYYY.mat';
       'CAFE', 'I:\Prj\npp_cafe_1M_1dg\grid_npp_nemo_1M_1dg_mat\npp_YYYY.mat'};
Tr = load(fullfile(fileparts(fileparts(OUT)), 'Extended_data_figures', 'ExtendedDataFig9_biome_sensitivity', ...
    'grid_monthly_trends_mld_npp_nemo.mat'), 'lat');
latG = double(Tr.lat(:));                         % row latitude of the (180,360) grid, exactly as interannual_biome.m
assert(numel(latG) == 180, 'unexpected latitude vector');
fprintf('grid latitude: first row %.1f, last row %.1f\n', latG(1), latG(end));
isN = latG >= 0;
bcell = cell(NB,1);
for b = 1:NB
    d = load(fullfile('I:\Prj\MLD_1M_NEMO\biomes_1M_1dg\mld\2003', ['mld_' biomes{b} '.mat']), 'idx');
    [ilon, ilat] = ind2sub([360 180], round(d.idx(:))); bcell{b} = sub2ind([180 360], ilat, ilon);
end
MLDm = nan(180,360,12,NY,'single');
for y = 1:NY, d = load(['I:\Prj\MLD_1M_NEMO\one_d03\mld\mld_' num2str(YEARS(y)) '.mat'], 'mldG1deg'); MLDm(:,:,:,y) = d.mldG1deg; end
MLDm = reshape(MLDm, 180*360, 12, NY);
pct = @(V) (V - mean(V, 2)) ./ mean(V, 2) * 100;     % rows = samples, cols = years
R = {};
for k = 1:3
    NPPm = nan(180,360,12,NY,'single');
    for y = 1:NY, d = load(strrep(src{k,2}, 'YYYY', num2str(YEARS(y))), 'nppG1dg'); NPPm(:,:,:,y) = d.nppG1dg; end
    NPPm = reshape(NPPm, 180*360, 12, NY);
    for s = 1:4
        Yall = []; Xall = []; Ball = []; nB = zeros(NB,1);
        for b = 1:NB
            c = bcell{b}; [ilat, ~] = ind2sub([180 360], c);
            for mi = 1:3
                for h = 1:2
                    if h == 1, cc = c(isN(ilat)); m = mmN(s,mi); else, cc = c(~isN(ilat)); m = mmS(s,mi); end
                    if isempty(cc), continue; end
                    Y = double(squeeze(NPPm(cc, m, :))); X = double(squeeze(MLDm(cc, m, :)));
                    if isvector(Y) && numel(cc) == 1, Y = Y(:)'; X = X(:)'; end
                    ok = all(isfinite(Y), 2) & all(isfinite(X), 2) & all(X > 0, 2) & all(Y >= 0, 2) & mean(Y, 2) > 0;
                    if ~any(ok), continue; end
                    Yd = detrend(pct(Y(ok,:))')'; Xd = detrend(pct(X(ok,:))')';
                    Yall = [Yall; Yd(:)]; Xall = [Xall; Xd(:)]; Ball = [Ball; repmat(b, numel(Yd), 1)]; %#ok<AGROW>
                    nB(b) = nB(b) + sum(ok);
                end
            end
        end
        T = table(Yall, Xall, categorical(Ball), 'VariableNames', {'lnNPP','lnMLD','biome'});
        out = nan(NB, 3); fx = nan(1,3);
        try
            mdl = fitlme(T, 'lnNPP ~ lnMLD + (1 + lnMLD | biome)'); ci = coefCI(mdl);
            b0 = mdl.Coefficients.Estimate(2); se0 = mdl.Coefficients.SE(2); fx = [b0 ci(2,:)];
            [~, ~, st] = randomEffects(mdl); isB = strcmp(st.Name, 'lnMLD');
            lev = str2double(st.Level(isB)); bb = st.Estimate(isB); sp = st.SEPred(isB);
            for j = 1:numel(lev), e = b0 + bb(j); hw = 1.96 * sqrt(se0^2 + sp(j)^2); out(lev(j),:) = [e e - hw e + hw]; end
        catch e, fprintf('  LME fail %s %s: %s\n', src{k,1}, seasName{s}, e.message);
        end
        R(end+1,:) = {src{k,1}, seasName{s}, 'All biomes', sum(nB), fx(1), fx(2), fx(3)}; %#ok<SAGROW>
        for b = 1:NB, R(end+1,:) = {src{k,1}, seasName{s}, biomes{b}, nB(b), out(b,1), out(b,2), out(b,3)}; end %#ok<SAGROW>
        fprintf('%s %s: %d samples (%d rows)\n', src{k,1}, seasName{s}, sum(nB), numel(Yall));
    end
end
A = cell2table(R, 'VariableNames', {'npp','season','biome','n','slope','ci_lo','ci_hi'});
A.slope(A.n == 0) = NaN; A.ci_lo(A.n == 0) = NaN; A.ci_hi(A.n == 0) = NaN;
writetable(A, fullfile(OUT, 'interannual_biome_allcells.csv'));
save(fullfile(OUT, 'interannual_biome_allcells.mat'), 'A', 'biomes', 'seasName');

F = load(fullfile(OUT, 'interannual_biome_pct.mat'), 'Tab'); F = F.Tab;
ordB = {'NA SPSS','NP SPSS','NA STSS','NP STSS','NA STPS','NP STPS','EP EQU','A EQU','WP EQU','IND STPS','SA STPS','SP STPS','SO STSS','SO SPSS'};
models = src(:,1)';
fprintf('\nALL CELLS vs FIGURE R4 SAMPLES: biome slope (n); flag - / + = CI excludes 0; biomes north to south\n');
for k = 1:3
    fprintf('-- %s\n', models{k});
    for un = [{'All biomes'}, ordB]
        fprintf('  %-10s', un{1});
        for s = 1:4
            a = A(strcmp(A.npp, models{k}) & strcmp(A.biome, un{1}) & strcmp(A.season, seasName{s}), :);
            f = F(strcmp(F.npp, models{k}) & strcmp(F.biome, un{1}) & strcmp(F.season, seasName{s}), :);
            fprintf(' | %5.2f%s(%5d) vs %5.2f%s(%4d)', a.slope, fl(a), a.n, f.slope, fl(f), f.n);
        end
        fprintf('\n');
    end
end
fprintf('\nSummary per model over 56 biome-seasons:\n');
for k = 1:3
    a = A(strcmp(A.npp, models{k}) & ~strcmp(A.biome, 'All biomes'), :);
    f = F(strcmp(F.npp, models{k}) & ~strcmp(F.biome, 'All biomes'), :);
    [~, ia, ib] = intersect(strcat(a.season, a.biome), strcat(f.season, f.biome));
    both = isfinite(a.slope(ia)) & isfinite(f.slope(ib));
    fprintf('  %s all cells: CI<0 %2d, CI>0 %2d, spans 0 %2d | Figure R4: CI<0 %2d, CI>0 %2d, spans 0 %2d, none %2d | same sign %d of %d\n', ...
        models{k}, sum(a.ci_hi < 0), sum(a.ci_lo > 0), sum(a.ci_lo <= 0 & a.ci_hi >= 0), ...
        sum(f.ci_hi < 0), sum(f.ci_lo > 0), sum(f.ci_lo <= 0 & f.ci_hi >= 0), sum(f.n == 0), ...
        sum(sign(a.slope(ia(both))) == sign(f.slope(ib(both)))), sum(both));
end
fprintf('\nEnd %s\n', datestr(now)); diary off

function s = fl(q)
    s = ' ';
    if isempty(q) || ~isfinite(q.slope), s = '?'; return; end
    if q.ci_hi < 0, s = '-'; elseif q.ci_lo > 0, s = '+'; end
end
