% Used for: carbon sensitivity, senMLD_carbon_v2.mat (Fig. R1)
% Original: Manuscript\MLD\manuscript\nature\NCC\R1\code for revision\sensitivty_carbon_mld\senMLD2_carbon_ot.m (data paths updated to I:\Prj; computation unchanged)
% Sensitivity of Carbon to MLD (all trends are included) significance level = 0.1
% Haipeng Zhao
% 3/20/2024
% 3/31/2024
% 1/28/2026
% 2/6/2026

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

senCarbon_NEMO = nan(6, 4, 3); % [lat_bands, seasons, statistics]
senCarbon_ECCO = nan(6, 4, 3);
senCarbon_HYCOM = nan(6, 4, 3);
senMLD_NEMO = nan(6, 4, 3);
senMLD_ECCO = nan(6, 4, 3);
senMLD_HYCOM = nan(6, 4, 3);

load CarbonTrends.mat

lat = flip( lat );
% carbon_trends_pct_nemo(carbon_trends_sig_nemo > 0.05) = nan;
% carbon_trends_pct_ecco(carbon_trends_sig_ecco > 0.05) = nan;
% carbon_trends_pct_hycom(carbon_trends_sig_hycom > 0.05) = nan;
% par_trends(par_trends_sig > 0.05) = nan;
% ws_trends(ws_trends_sig > 0.05) = nan;
% sla_trends(sla_trends_sig > 0.05) = nan;

for i = 1 : 6
    latIdx = lat < latBand(i, 1) & lat > latBand(i, 2);
    for j = 1 : 4
        %% NEMO
        c1_carbon = nan( 120, 3 ); c1_carbon_sig = nan( 120, 3 ); % Carbon_NEMO
        c1_carbon(lat > 0, :) = carbon_trends_pct_nemo(lat > 0, mmN(j, :)); c1_carbon_sig(lat > 0, :) = carbon_trends_sig_nemo(lat > 0, mmN(j, :));
        c1_carbon(lat < 0, :) = carbon_trends_pct_nemo( lat < 0, mmS(j, :)); c1_carbon_sig(lat < 0, :) = carbon_trends_sig_nemo( lat < 0, mmS(j, :));
        c1_carbon = c1_carbon(latIdx, :); c1_carbon_sig = c1_carbon_sig(latIdx, :);
        c1_carbon = c1_carbon(:); c1_carbon_sig = c1_carbon_sig(:);
        idx_carbon_sig = ~isnan( c1_carbon );
        c1_carbon = c1_carbon( idx_carbon_sig );
        c1_carbon_sig = c1_carbon_sig( idx_carbon_sig );

        c1_mld = nan( 120, 3 ); c1_mld_sig = nan( 120, 3 ); % MLD_NEMO
        c1_mld(lat > 0, :) = mld_trends_pct_nemo(lat > 0, mmN(j, :)); c1_mld_sig(lat > 0, :) = mld_trends_sig_nemo(lat > 0, mmN(j, :));
        c1_mld(lat < 0, :) = mld_trends_pct_nemo( lat < 0, mmS(j, :)); c1_mld_sig(lat < 0, :) = mld_trends_sig_nemo( lat < 0, mmS(j, :));
        c1_mld = c1_mld(latIdx, :); c1_mld_sig = c1_mld_sig(latIdx, :);
        c1_mld = c1_mld(:); c1_mld_sig = c1_mld_sig(:);
        c1_mld = c1_mld( idx_carbon_sig );
        c1_mld_sig = c1_mld_sig( idx_carbon_sig );

        %% ECCO
        c2_carbon = nan( 120, 3 ); c2_carbon_sig = nan( 120, 3 ); % Carbon_ECCO
        c2_carbon(lat > 0, :) = carbon_trends_pct_ecco(lat > 0, mmN(j, :)); c2_carbon_sig(lat > 0, :) = carbon_trends_sig_ecco(lat > 0, mmN(j, :));
        c2_carbon(lat < 0, :) = carbon_trends_pct_ecco( lat < 0, mmS(j, :)); c2_carbon_sig(lat < 0, :) = carbon_trends_sig_ecco( lat < 0, mmS(j, :));
        c2_carbon = c2_carbon(latIdx, :); c2_carbon_sig = c2_carbon_sig(latIdx, :);
        c2_carbon = c2_carbon(:); c2_carbon_sig = c2_carbon_sig(:);
        idx_carbon_sig = ~isnan( c2_carbon );
        c2_carbon = c2_carbon( idx_carbon_sig );
        c2_carbon_sig = c2_carbon_sig( idx_carbon_sig );

        c2_mld = nan( 120, 3 ); c2_mld_sig = nan( 120, 3 ); % MLD_ECCO
        c2_mld(lat > 0, :) = mld_trends_pct_ecco(lat > 0, mmN(j, :)); c2_mld_sig(lat > 0, :) = mld_trends_sig_ecco(lat > 0, mmN(j, :));
        c2_mld(lat < 0, :) = mld_trends_pct_ecco( lat < 0, mmS(j, :)); c2_mld_sig(lat < 0, :) = mld_trends_sig_ecco( lat < 0, mmS(j, :));
        c2_mld = c2_mld(latIdx, :); c2_mld_sig = c2_mld_sig(latIdx, :);
        c2_mld = c2_mld(:); c2_mld_sig = c2_mld_sig(:);
        c2_mld = c2_mld( idx_carbon_sig );
        c2_mld_sig = c2_mld_sig( idx_carbon_sig );

        %% HYCOM
        c3_carbon = nan( 120, 3 ); c3_carbon_sig = nan( 120, 3 ); % Carbon_HYCOM
        c3_carbon(lat > 0, :) = carbon_trends_pct_hycom(lat > 0, mmN(j, :)); c3_carbon_sig(lat > 0, :) = carbon_trends_sig_hycom(lat > 0, mmN(j, :));
        c3_carbon(lat < 0, :) = carbon_trends_pct_hycom( lat < 0, mmS(j, :)); c3_carbon_sig(lat < 0, :) = carbon_trends_sig_hycom( lat < 0, mmS(j, :));
        c3_carbon = c3_carbon(latIdx, :); c3_carbon_sig = c3_carbon_sig(latIdx, :);
        c3_carbon = c3_carbon(:); c3_carbon_sig = c3_carbon_sig(:);
        idx_carbon_sig = (~isnan( c3_carbon ));
        c3_carbon = c3_carbon( idx_carbon_sig );
        c3_carbon_sig = c3_carbon_sig( idx_carbon_sig );

        c3_mld = nan( 120, 3 ); c3_mld_sig = nan( 120, 3 ); % MLD_HYCOM
        c3_mld(lat > 0, :) = mld_trends_pct_hycom(lat > 0, mmN(j, :)); c3_mld_sig(lat > 0, :) = mld_trends_sig_hycom(lat > 0, mmN(j, :));
        c3_mld(lat < 0, :) = mld_trends_pct_hycom( lat < 0, mmS(j, :)); c3_mld_sig(lat < 0, :) = mld_trends_sig_hycom( lat < 0, mmS(j, :));
        c3_mld = c3_mld(latIdx, :); c3_mld_sig = c3_mld_sig(latIdx, :);
        c3_mld = c3_mld(:); c3_mld_sig = c3_mld_sig(:);
        c3_mld = c3_mld( idx_carbon_sig );
        c3_mld_sig = c3_mld_sig( idx_carbon_sig );

        %%
        sig_level = 0.1;
        idx_sig = (c1_carbon_sig < sig_level) & (c1_mld_sig < sig_level);
        sen_NEMO = c1_carbon(idx_sig) ./ c1_mld(idx_sig);
        if ~isempty( sen_NEMO )
           senMLD_NEMO(i, j, 1) = nanmedian( sen_NEMO(:, 1) );
           senMLD_NEMO(i, j, 2) = nanstd( sen_NEMO(:, 1) );
           senMLD_NEMO(i, j, 3) = sum( ~isnan( sen_NEMO(:, 1) ) );
        end
        c1_mld(~idx_sig) = 0;
        senPt_nemo{i , j} = [c1_carbon, c1_mld];

        idx_sig = (c2_carbon_sig < sig_level) & (c2_mld_sig < sig_level);
        sen_ECCO = c2_carbon(idx_sig) ./ c2_mld(idx_sig);
        if ~isempty( sen_ECCO )
           senMLD_ECCO(i, j, 1) = nanmedian( sen_ECCO(:, 1) );
           senMLD_ECCO(i, j, 2) = nanstd( sen_ECCO(:, 1) );
           senMLD_ECCO(i, j, 3) = sum( ~isnan( sen_ECCO(:, 1) ) );
        end
        c2_mld(~idx_sig) = 0;
        senPt_ecco{i , j} = [c2_carbon, c2_mld];

       idx_sig = (c3_carbon_sig < sig_level) & (c3_mld_sig < sig_level);
       sen_HYCOM = c3_carbon(idx_sig) ./ c3_mld(idx_sig);
       if ~isempty( sen_HYCOM )
           senMLD_HYCOM(i, j, 1) = nanmedian( sen_HYCOM(:, 1) );
           senMLD_HYCOM(i, j, 2) = nanstd( sen_HYCOM(:, 1) );
           senMLD_HYCOM(i, j, 3) = sum( ~isnan( sen_HYCOM(:, 1) ) );
       end
        c3_mld(~idx_sig) = 0;
        senPt_hycom{i , j} = [c3_carbon, c3_mld];

    end
end

save( 'senMLD_carbon_v2.mat', 'senMLD_NEMO', 'senMLD_ECCO', 'senMLD_HYCOM', 'senPt_nemo', 'senPt_ecco', 'senPt_hycom' )
