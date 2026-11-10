% Used for: driver for Fig. R1 and the sample-count comparison
% Original: Manuscript\MLD\manuscript\nature\NCC\R2\code for revision\fig4_p005\run_fig4_p005.m (data paths updated to I:\Prj; computation unchanged)
% Reviewer 1 (R2 round): Fig. 4 at p<0.05. Runs the copies, checks that the p<0.1 copy
% reproduces the published Code\senMLD_v2.mat, then plots p<0.05.
cd(fileparts(mfilename('fullpath')))
diary('run_fig4_p005_log.txt'); diary on

senMLD2_p01; senMLD2_p005;
A = load('senMLD_v2_p01.mat');
B = load(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), '4_trends_and_sensitivity', 'senMLD_v2.mat'));
f = {'senMLD_NEMO','senMLD_ECCO','senMLD_HYCOM'};
for k = 1:3
    d = abs(A.(f{k}) - B.(f{k})); d(isnan(A.(f{k})) & isnan(B.(f{k}))) = 0;
    fprintf('p<0.1 copy vs published %s: max diff %g, NaN pattern equal %d\n', f{k}, max(d(:)), isequal(isnan(A.(f{k})), isnan(B.(f{k}))));
end

C = load('senMLD_v2_p005.mat');
fprintf('\nBand-season cells (24 per product): neg / pos / empty, and neg / pos with n>=3\n');
for k = 1:3
    for P = {A, C}
        m = P{1}.(f{k})(:,:,1); n = P{1}.(f{k})(:,:,3);
        fprintf('%-13s neg %2d pos %2d empty %2d | n>=3: neg %2d pos %2d\n', f{k}, sum(m(:)<0), sum(m(:)>0), sum(isnan(m(:))), sum(m(:)<0 & n(:)>=3), sum(m(:)>0 & n(:)>=3));
    end
end
fprintf('(first row of each product = p<0.1, second = p<0.05)\n');
fprintf('\nNEMO median [n] per band (rows 60-40N ... 40-60S) x season (Sp Su Au Wi):\n');
for P = {A, C}
    m = P{1}.senMLD_NEMO(:,:,1); n = P{1}.senMLD_NEMO(:,:,3);
    for i = 1:6, fprintf('%7.2f[%3d]', [m(i,:); n(i,:)]); fprintf('\n'); end
    fprintf('\n');
end

f4_plot_sensitivity_p005;
diary off
