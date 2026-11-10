% Used for: climatological mixed-layer PAR and MLD, par_mld_mldonly.mat (Fig. S11)
% Original: Manuscript\MLD\manuscript\nature\NCC\R1\code for revision\mld_npp_relationship_grid\compute_par_mld_mldonly.m (data paths updated to I:\Prj; computation unchanged)
% compute_par_mld_mldonly.m
%
% Rebuilds CbPM's mixed-layer light climate par_mld from its three inputs
% (surface PAR, Kd490, MLD) so that the MLD-driven part of the light change
% can be isolated:
%
%   par_mld = [ int_lambda Eo(lambda)*0.975*exp(-K(lambda)*MLD/2) ] / daylength_h
%   Eo(lambda) = PAR * parFraction(lambda),  K(lambda) = austinPetzold_1986(lambda, Kd490)
%
% Three fields are produced, all with PAR and Kd490 held at their 2003-2022
% monthly climatology so that ONLY MLD varies:
%   par_mld_clim : MLD at its 2003-2022 monthly climatology
%   par_mld_2003 : MLD at its 2003 Theil-Sen trend endpoint
%   par_mld_2022 : MLD at its 2022 Theil-Sen trend endpoint
% so (par_mld_2022 - par_mld_2003) is the light change attributable to the
% MLD trend alone.
%
% A fourth field, par_mld_allvary, uses each year's actual PAR/Kd490/MLD and
% is kept only to quantify how much of the total light change is MLD-driven.
%
% ORDER OF OPERATIONS: the formula is nonlinear in MLD and Kd, so it is
% applied per pixel per calendar month and only then averaged. A previous
% attempt in this project averaged inputs first and was invalid; see the
% 2026-08-23 entry in revision_conversation_log.md.
%
% Two validations run before anything is produced and the script aborts if
% either fails:
%   V1 reproduce the stored par_mld for 2003 from that year's own inputs
%   V2 vectorised Theil-Sen endpoints against ktaub on a random sample
%
% Output: par_mld_mldonly.mat (this folder)
%
% Haipeng Zhao
% 2026-09-10

clear; clc;
OUT_DIR = fileparts(mfilename('fullpath'));
CBPM_DIR = fullfile(fileparts(fileparts(OUT_DIR)), 'helper_functions');
addpath(CBPM_DIR);
addpath('D:\Users\zhp19\OneDrive - Indiana State University\Research\MATLAB toolbox\ktaub');

YEARS = 2003:2022; NY = numel(YEARS);
LAMBDA      = [400, 412, 443, 490, 510, 555, 625, 670, 700];
parFraction = [0.0029, 0.0032, 0.0035, 0.0037, 0.0037, 0.0036, 0.0032, 0.0030, 0.0024];

%% ---------- load inputs ----------
fprintf('Loading inputs ...\n');
fp = dir('I:\Prj\npp_cbpm_1M_1dg\grid_par\*.mat');
fk = dir('I:\Prj\npp_cbpm_1M_1dg\grid_kd490\*.mat');
fm = dir('I:\Prj\MLD_1M_NEMO\one_d03\mld\*.mat');

PAR = nan(180,360,12,NY,'single');
KD  = nan(180,360,12,NY,'single');
MLD = nan(180,360,12,NY,'single');
for y = 1:NY
    s = num2str(YEARS(y));
    a = load(fullfile(fp(find(contains({fp.name},s),1)).folder, fp(find(contains({fp.name},s),1)).name));
    b = load(fullfile(fk(find(contains({fk.name},s),1)).folder, fk(find(contains({fk.name},s),1)).name));
    c = load(fullfile(fm(find(contains({fm.name},s),1)).folder, fm(find(contains({fm.name},s),1)).name));
    PAR(:,:,:,y) = a.par;
    fn = fieldnames(b); KD(:,:,:,y)  = b.(fn{find(~ismember(fn,{'lat','lon'}),1)});
    fn = fieldnames(c); MLD(:,:,:,y) = c.(fn{find(contains(fn,'mld','IgnoreCase',true),1)});
end
ref = load('I:\Prj\npp_cbpm_1M_1dg_zno3mMU\grid_npp_nemo\npp_2003.mat','par_mld','lat');
lat = double(ref.lat(:));

%% ---------- daylength, exactly as CbPM computes it ----------
DLh = nan(180,12);
for m = 1:12
    nd = eomday(2003,m); nDays = yeardays(2003); acc = zeros(180,nd);
    for i = 1:nd
        yd = day(datetime(2003,m,i),'dayofyear');
        decl = 23.5*cos(2*pi*(yd-172)/nDays)*pi/180;
        v = -tan(lat*pi/180).*tan(decl);
        acc(:,i) = acos(min(max(v,-1),1))/pi;
    end
    DLh(:,m) = mean(acc,2)*24;   % decimal hours
end

%% ---------- V1: reproduce stored par_mld for 2003 ----------
fprintf('V1 reproducing stored par_mld (2003) ...\n');
chk = par_mld_from(PAR(:,:,:,1), KD(:,:,:,1), MLD(:,:,:,1), LAMBDA, parFraction, DLh);
r = double(ref.par_mld); c = double(chk);
v = isfinite(r) & isfinite(c) & r > 0;
rel = abs(c(v)-r(v))./r(v);
fprintf('   n=%d  median rel diff=%.3e  p99=%.3e  max=%.3e\n', sum(v), median(rel), prctile(rel,99), max(rel));
if median(rel) > 1e-3
    error('V1 FAILED: reimplementation does not reproduce stored par_mld.');
end
fprintf('   V1 passed.\n');

%% ---------- Theil-Sen endpoints per pixel per month (vectorised) ----------
fprintf('Computing Theil-Sen MLD endpoints ...\n');
x = (1:NY)';
[ia,ib] = find(triu(ones(NY),1));       % 190 pairs
MLD2003 = nan(180,360,12,'single');
MLD2022 = nan(180,360,12,'single');
for m = 1:12
    Y = reshape(double(MLD(:,:,m,:)), [], NY)';        % NY x Npix
    ok = all(isfinite(Y),1);
    Yv = Y(:,ok);
    slopes = (Yv(ib,:)-Yv(ia,:)) ./ (x(ib)-x(ia));
    sen = median(slopes,1);
    % ktaub anchors the Sen line at datain(round(n/2),1), i.e. x(10) for
    % n=20, NOT at median(x)=10.5. Matching that exactly (ktaub.m line 391).
    inter = median(Yv,1) - sen.*x(round(NY/2));
    a = nan(1,size(Y,2)); b = nan(1,size(Y,2));
    a(ok) = inter + sen*x(1);
    b(ok) = inter + sen*x(end);
    MLD2003(:,:,m) = single(reshape(a,180,360));
    MLD2022(:,:,m) = single(reshape(b,180,360));
end

%% ---------- V2: vectorised Sen vs ktaub on a sample ----------
fprintf('V2 checking Sen endpoints against ktaub ...\n');
rng(0); nchk = 40; d2003 = nan(nchk,1); d2022 = nan(nchk,1); got = 0; tries = 0;
while got < nchk && tries < 4000
    tries = tries + 1;
    i = randi(180); j = randi(360); m = randi(12);
    y = squeeze(double(MLD(i,j,m,:)));
    if any(~isfinite(y)) || range(y) == 0, continue; end
    [~,~,~,~,~,~,~,~,~,pos,~] = ktaub([x, y], 0.9999, 1);
    if any(~isfinite(pos(:))), continue; end
    got = got + 1;
    d2003(got) = abs(pos(1,2)   - double(MLD2003(i,j,m)));
    d2022(got) = abs(pos(end,2) - double(MLD2022(i,j,m)));
end
fprintf('   sampled %d pixels: max |diff| 2003 = %.3e m, 2022 = %.3e m\n', got, max(d2003), max(d2022));
% MLD is stored as single, so ~1e-5 m on a ~50 m value is rounding, not error.
if max([d2003;d2022]) > 1e-3
    error('V2 FAILED: vectorised Sen endpoints disagree with ktaub.');
end
fprintf('   V2 passed.\n');

%% ---------- climatologies and the three counterfactual fields ----------
PAR_clim = mean(PAR,4,'omitnan');
KD_clim  = mean(KD ,4,'omitnan');
MLD_clim = mean(MLD,4,'omitnan');

fprintf('Building counterfactual par_mld fields ...\n');
par_mld_clim = par_mld_from(PAR_clim, KD_clim, MLD_clim, LAMBDA, parFraction, DLh);
par_mld_2003 = par_mld_from(PAR_clim, KD_clim, MLD2003 , LAMBDA, parFraction, DLh);
par_mld_2022 = par_mld_from(PAR_clim, KD_clim, MLD2022 , LAMBDA, parFraction, DLh);

% all-varying reference, for attribution only
par_mld_allvary = nan(180,360,12,NY,'single');
for y = 1:NY
    par_mld_allvary(:,:,:,y) = par_mld_from(PAR(:,:,:,y), KD(:,:,:,y), MLD(:,:,:,y), LAMBDA, parFraction, DLh);
end

save(fullfile(OUT_DIR,'par_mld_mldonly.mat'), ...
    'par_mld_clim','par_mld_2003','par_mld_2022','par_mld_allvary', ...
    'MLD_clim','MLD2003','MLD2022','PAR_clim','KD_clim','lat','YEARS','DLh','-v7.3');
fprintf('Saved par_mld_mldonly.mat\n');

d = double(par_mld_2022) - double(par_mld_2003);
fprintf('MLD-only par_mld change (mol/m2/h): median %+.4f, p5 %+.4f, p95 %+.4f\n', ...
    median(d(isfinite(d))), prctile(d(isfinite(d)),5), prctile(d(isfinite(d)),95));

%% ---------- local functions ----------
function out = par_mld_from(parF, kdF, mldF, LAMBDA, parFraction, DLh)
out = nan(size(parF),'single');
for mm = 1:12
    p = double(parF(:,:,mm)); k = double(kdF(:,:,mm)); d = double(mldF(:,:,mm));
    acc2 = zeros(180,360); Ezprev = [];
    for q = 1:numel(LAMBDA)
        Kq  = austinPetzold_1986(LAMBDA(q), k);
        Ezq = p .* parFraction(q) .* 0.975 .* exp(-Kq .* d / 2);
        if q > 1
            acc2 = acc2 + (LAMBDA(q)-LAMBDA(q-1)) .* (Ezq + Ezprev)/2;
        end
        Ezprev = Ezq;
    end
    out(:,:,mm) = single(acc2 ./ repmat(DLh(:,mm),1,360));
end
end
