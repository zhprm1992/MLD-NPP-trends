% Used for: trends per latitude row and month for ensemble NPP and MLD-NEMO (SI Fig. S3)
% Original: Manuscript\MLD\manuscript\nature\NCC\R1\code for revision\npp_ensemble\npp_mld_trends_ensemble_allProducts_ot.m (data paths updated to I:\Prj; computation unchanged)
% Compute Mann-Kendall trends (% yr-1) of NPP (ensemble: VGPM+CAFE+CbPM)
% and MLD, by latitude and month, for all 3 MLD products (NEMO/ECCO/HYCOM)
% -- generalizes npp_mld_trends_ensemble_ot.m (NEMO-only) to loop over
% products and save into one combined file with per-product suffixed
% variables, matching CAFE_NPPtrends.mat's structure exactly
% (npp_trends_pct_nemo, npp_trends_pct_ecco, npp_trends_pct_hycom, ...).
%
% Haipeng Zhao
% 2026-07-15

clear
clc

products = {'nemo','ecco','hycom'};
mld_files = { ...
    'I:\Prj\MLD_1M_NEMO\mldLat_monthly_nemo.mat', ...
    'I:\Prj\MLD_ECCO\mldLat_monthly_ecco.mat', ...
    'I:\Prj\MLD_HYCOM\mldLat_monthly_hycom.mat'};

results = struct();

for pi = 1:numel(products)
    prod = products{pi};
    fprintf('=== Product: %s ===\n', upper(prod));

    d = load( fullfile('I:\Prj\npp_ensemble_1M_1dg', ['nppLat_monthly_' prod '.mat']) );
    idx = d.lat < 60 & d.lat > -60;
    lat = d.lat(idx);
    nppLat = d.nppLat(idx, :);
    mm_npp = d.mm;

    d2 = load( mld_files{pi} );
    idx2 = d2.lat < 60 & d2.lat > -60;
    mldLat = d2.mldLat(idx2, :);
    mm_mld = d2.mm;

    npp_trends_pct = nan( 12, length( lat ) );
    npp_trends_sig = nan( 12, length( lat ) );
    mld_trends_pct = nan( 12, length( lat ) );
    mld_trends_sig = nan( 12, length( lat ) );

    %% NPP trends
    for i = 1 : length( lat )
        for j = 1 : 12
            npp_lat = nppLat(i, mm_npp == j);
            if sum( isnan( npp_lat ) ) > 0
                continue
            end
            period_length = size( npp_lat, 2 );
            datain = [(1 : period_length)', npp_lat'];
            [~, ~, r] = regress( datain(:, 2), [(ones( length( datain(:, 1) ), 1 )), datain(:, 1)]  );
            dwtest( r, datain(:, 1) );
            [~, ~, ~, sig, ~, ~, ~, sen, ~, plotofslope, ~] = ktaub(datain, 0.9999, 1);
            p_value = sig;
            if ~isnan(plotofslope)
                if plotofslope(1, 2) == 0
                    npp_trends_pct(j, i) = 0;
                else
                    npp_trends_pct(j, i) = single( (plotofslope(2, 2) - plotofslope(1, 2)) / plotofslope(1, 2) * 100 );
                end
                npp_trends_sig(j, i) = p_value;
            end
        end
    end

    %% MLD trends
    for i = 1 : length( lat )
        for j = 1 : 12
            mld_lat = mldLat(i, mm_mld == j);
            if sum( isnan( mld_lat ) ) > 0
                continue
            end
            period_length = size( mld_lat, 2 );
            datain = [(1 : period_length)', mld_lat'];
            [~, ~, r] = regress( datain(:, 2), [(ones( length( datain(:, 1) ), 1 )), datain(:, 1)]  );
            dwtest( r, datain(:, 1) );
            [~, ~, ~, sig, ~, ~, ~, sen, ~, plotofslope, ~] = ktaub(datain, 0.9999, 1);
            p_value = sig;
            if ~isnan(plotofslope)
                if plotofslope(1, 2) == 0
                    mld_trends_pct(j, i) = 0;
                else
                    mld_trends_pct(j, i) = single( (plotofslope(2, 2) - plotofslope(1, 2)) / plotofslope(1, 2) * 100 );
                end
                mld_trends_sig(j, i) = p_value;
            end
        end
    end

    % Rotate to [lat x month]
    results.(['npp_trends_pct_' prod]) = rot90( npp_trends_pct );
    results.(['npp_trends_sig_' prod]) = rot90( npp_trends_sig );
    results.(['mld_trends_pct_' prod]) = rot90( mld_trends_pct );
    results.(['mld_trends_sig_' prod]) = rot90( mld_trends_sig );
    results.lat = lat;
end

npp_trends_pct_nemo = results.npp_trends_pct_nemo; npp_trends_sig_nemo = results.npp_trends_sig_nemo;
mld_trends_pct_nemo = results.mld_trends_pct_nemo; mld_trends_sig_nemo = results.mld_trends_sig_nemo;
npp_trends_pct_ecco = results.npp_trends_pct_ecco; npp_trends_sig_ecco = results.npp_trends_sig_ecco;
mld_trends_pct_ecco = results.mld_trends_pct_ecco; mld_trends_sig_ecco = results.mld_trends_sig_ecco;
npp_trends_pct_hycom = results.npp_trends_pct_hycom; npp_trends_sig_hycom = results.npp_trends_sig_hycom;
mld_trends_pct_hycom = results.mld_trends_pct_hycom; mld_trends_sig_hycom = results.mld_trends_sig_hycom;
lat = results.lat;

save( 'I:\Prj\npp_ensemble_1M_1dg\ENSEMBLE_NPPtrends.mat', ...
      'npp_trends_pct_nemo', 'mld_trends_pct_nemo', 'npp_trends_sig_nemo', 'mld_trends_sig_nemo', ...
      'npp_trends_pct_ecco', 'mld_trends_pct_ecco', 'npp_trends_sig_ecco', 'mld_trends_sig_ecco', ...
      'npp_trends_pct_hycom', 'mld_trends_pct_hycom', 'npp_trends_sig_hycom', 'mld_trends_sig_hycom', ...
      'lat' )
disp( 'Save successfully: ENSEMBLE_NPPtrends (all 3 products)' )
