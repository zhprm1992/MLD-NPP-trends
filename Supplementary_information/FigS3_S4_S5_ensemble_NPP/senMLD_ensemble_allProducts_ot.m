% Used for: ensemble sensitivity per latitude band and season, senMLD_ensemble.mat (SI Fig. S3)
% Original: Manuscript\MLD\manuscript\nature\NCC\R1\code for revision\npp_ensemble\senMLD_ensemble_allProducts_ot.m (data paths updated to I:\Prj; computation unchanged)
% Sensitivity of NPP (ensemble: VGPM+CAFE+CbPM) to MLD, all 3 products
% (NEMO/ECCO/HYCOM), significance level = 0.1 -- generalizes
% senMLD_ensemble_ot.m (NEMO-only) to loop over products and save with
% per-product suffixed variables, matching npp_cafe\senMLD2_ot.m's
% structure (senMLD_NEMO/ECCO/HYCOM, senPt_nemo/ecco/hycom).
%
% Haipeng Zhao
% 2026-07-15

clear
clc

mmN = [3 4 5;
                    6 7 8;
                    9 10 11;
                    12 1 2; ];  % [Spring, Summer, Autumn, Winter]

mmS = [9 10 11;
                 12 1 2;
                   3 4 5;
                   6 7 8];

latBand = [60 40;
                      40 20;
                      20 0;
                      0 -20;
                      -20 -40;
                      -40 -60];

senMLD_NEMO = nan(6, 4, 3);
senMLD_ECCO = nan(6, 4, 3);
senMLD_HYCOM = nan(6, 4, 3);

load I:\Prj\npp_ensemble_1M_1dg\ENSEMBLE_NPPtrends.mat

lat = flip( lat );

for i = 1 : 6
    latIdx = lat < latBand(i, 1) & lat > latBand(i, 2);
    for j = 1 : 4
        %% NEMO
        c1_npp = nan( 120, 3 ); c1_npp_sig = nan( 120, 3 );
        c1_npp(lat > 0, :) = npp_trends_pct_nemo(lat > 0, mmN(j, :)); c1_npp_sig(lat > 0, :) = npp_trends_sig_nemo(lat > 0, mmN(j, :));
        c1_npp(lat < 0, :) = npp_trends_pct_nemo( lat < 0, mmS(j, :)); c1_npp_sig(lat < 0, :) = npp_trends_sig_nemo( lat < 0, mmS(j, :));
        c1_npp = c1_npp(latIdx, :); c1_npp_sig = c1_npp_sig(latIdx, :);
        c1_npp = c1_npp(:); c1_npp_sig = c1_npp_sig(:);
        idx_NPP_sig = ~isnan( c1_npp );
        c1_npp = c1_npp( idx_NPP_sig );
        c1_npp_sig = c1_npp_sig( idx_NPP_sig );

        c1_mld = nan( 120, 3 ); c1_mld_sig = nan( 120, 3 );
        c1_mld(lat > 0, :) = mld_trends_pct_nemo(lat > 0, mmN(j, :)); c1_mld_sig(lat > 0, :) = mld_trends_sig_nemo(lat > 0, mmN(j, :));
        c1_mld(lat < 0, :) = mld_trends_pct_nemo( lat < 0, mmS(j, :)); c1_mld_sig(lat < 0, :) = mld_trends_sig_nemo( lat < 0, mmS(j, :));
        c1_mld = c1_mld(latIdx, :); c1_mld_sig = c1_mld_sig(latIdx, :);
        c1_mld = c1_mld(:); c1_mld_sig = c1_mld_sig(:);
        c1_mld = c1_mld( idx_NPP_sig );
        c1_mld_sig = c1_mld_sig( idx_NPP_sig );

        %% ECCO
        c2_npp = nan( 120, 3 ); c2_npp_sig = nan( 120, 3 );
        c2_npp(lat > 0, :) = npp_trends_pct_ecco(lat > 0, mmN(j, :)); c2_npp_sig(lat > 0, :) = npp_trends_sig_ecco(lat > 0, mmN(j, :));
        c2_npp(lat < 0, :) = npp_trends_pct_ecco( lat < 0, mmS(j, :)); c2_npp_sig(lat < 0, :) = npp_trends_sig_ecco( lat < 0, mmS(j, :));
        c2_npp = c2_npp(latIdx, :); c2_npp_sig = c2_npp_sig(latIdx, :);
        c2_npp = c2_npp(:); c2_npp_sig = c2_npp_sig(:);
        idx_NPP_sig2 = ~isnan( c2_npp );
        c2_npp = c2_npp( idx_NPP_sig2 );
        c2_npp_sig = c2_npp_sig( idx_NPP_sig2 );

        c2_mld = nan( 120, 3 ); c2_mld_sig = nan( 120, 3 );
        c2_mld(lat > 0, :) = mld_trends_pct_ecco(lat > 0, mmN(j, :)); c2_mld_sig(lat > 0, :) = mld_trends_sig_ecco(lat > 0, mmN(j, :));
        c2_mld(lat < 0, :) = mld_trends_pct_ecco( lat < 0, mmS(j, :)); c2_mld_sig(lat < 0, :) = mld_trends_sig_ecco( lat < 0, mmS(j, :));
        c2_mld = c2_mld(latIdx, :); c2_mld_sig = c2_mld_sig(latIdx, :);
        c2_mld = c2_mld(:); c2_mld_sig = c2_mld_sig(:);
        c2_mld = c2_mld( idx_NPP_sig2 );
        c2_mld_sig = c2_mld_sig( idx_NPP_sig2 );

        %% HYCOM
        c3_npp = nan( 120, 3 ); c3_npp_sig = nan( 120, 3 );
        c3_npp(lat > 0, :) = npp_trends_pct_hycom(lat > 0, mmN(j, :)); c3_npp_sig(lat > 0, :) = npp_trends_sig_hycom(lat > 0, mmN(j, :));
        c3_npp(lat < 0, :) = npp_trends_pct_hycom( lat < 0, mmS(j, :)); c3_npp_sig(lat < 0, :) = npp_trends_sig_hycom( lat < 0, mmS(j, :));
        c3_npp = c3_npp(latIdx, :); c3_npp_sig = c3_npp_sig(latIdx, :);
        c3_npp = c3_npp(:); c3_npp_sig = c3_npp_sig(:);
        idx_NPP_sig3 = ~isnan( c3_npp );
        c3_npp = c3_npp( idx_NPP_sig3 );
        c3_npp_sig = c3_npp_sig( idx_NPP_sig3 );

        c3_mld = nan( 120, 3 ); c3_mld_sig = nan( 120, 3 );
        c3_mld(lat > 0, :) = mld_trends_pct_hycom(lat > 0, mmN(j, :)); c3_mld_sig(lat > 0, :) = mld_trends_sig_hycom(lat > 0, mmN(j, :));
        c3_mld(lat < 0, :) = mld_trends_pct_hycom( lat < 0, mmS(j, :)); c3_mld_sig(lat < 0, :) = mld_trends_sig_hycom( lat < 0, mmS(j, :));
        c3_mld = c3_mld(latIdx, :); c3_mld_sig = c3_mld_sig(latIdx, :);
        c3_mld = c3_mld(:); c3_mld_sig = c3_mld_sig(:);
        c3_mld = c3_mld( idx_NPP_sig3 );
        c3_mld_sig = c3_mld_sig( idx_NPP_sig3 );

        %%
        sig_level = 0.1;
        idx_sig = (c1_npp_sig < sig_level) & (c1_mld_sig < sig_level);
        sen_NEMO = c1_npp(idx_sig) ./ c1_mld(idx_sig);
        if ~isempty( sen_NEMO )
           senMLD_NEMO(i, j, 1) = nanmedian( sen_NEMO(:, 1) );
           senMLD_NEMO(i, j, 2) = nanstd( sen_NEMO(:, 1) );
           senMLD_NEMO(i, j, 3) = sum( ~isnan( sen_NEMO(:, 1) ) );
        end
        c1_mld(~idx_sig) = 0;
        senPt_nemo{i , j} = [c1_npp, c1_mld];

        idx_sig = (c2_npp_sig < sig_level) & (c2_mld_sig < sig_level);
        sen_ECCO = c2_npp(idx_sig) ./ c2_mld(idx_sig);
        if ~isempty( sen_ECCO )
           senMLD_ECCO(i, j, 1) = nanmedian( sen_ECCO(:, 1) );
           senMLD_ECCO(i, j, 2) = nanstd( sen_ECCO(:, 1) );
           senMLD_ECCO(i, j, 3) = sum( ~isnan( sen_ECCO(:, 1) ) );
        end
        c2_mld(~idx_sig) = 0;
        senPt_ecco{i , j} = [c2_npp, c2_mld];

       idx_sig = (c3_npp_sig < sig_level) & (c3_mld_sig < sig_level);
       sen_HYCOM = c3_npp(idx_sig) ./ c3_mld(idx_sig);
       if ~isempty( sen_HYCOM )
           senMLD_HYCOM(i, j, 1) = nanmedian( sen_HYCOM(:, 1) );
           senMLD_HYCOM(i, j, 2) = nanstd( sen_HYCOM(:, 1) );
           senMLD_HYCOM(i, j, 3) = sum( ~isnan( sen_HYCOM(:, 1) ) );
       end
        c3_mld(~idx_sig) = 0;
        senPt_hycom{i , j} = [c3_npp, c3_mld];
    end
end

save( 'I:\Prj\npp_ensemble_1M_1dg\senMLD_ensemble.mat', 'senMLD_NEMO', 'senMLD_ECCO', 'senMLD_HYCOM', 'senPt_nemo', 'senPt_ecco', 'senPt_hycom' )
disp( 'Save successfully: senMLD_ensemble (all 3 products)' )
