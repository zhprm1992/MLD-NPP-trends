% Used for: Fig. 4 sensitivity with p < 0.05 (Fig. R1)
% Original: Manuscript\MLD\manuscript\nature\NCC\R2\code for revision\fig4_p005\senMLD2_p005.m (data paths updated to I:\Prj; computation unchanged)
% R2 copy of Code\senMLD2_ot.m, only sig_level and file paths changed. significance level = 0.05
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

senNPP_NEMO = nan(6, 4, 3); % [lat_bands, seasons, statistics]
senNPP_ECCO = nan(6, 4, 3);
senNPP_HYCOM = nan(6, 4, 3);
senMLD_NEMO = nan(6, 4, 3);
senMLD_ECCO = nan(6, 4, 3);
senMLD_HYCOM = nan(6, 4, 3);

load(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), '4_trends_and_sensitivity', 'NPPtrends.mat'))

lat = flip( lat );
% npp_trends_pct_nemo(npp_trends_sig_nemo > 0.05) = nan;
% npp_trends_pct_ecco(npp_trends_sig_ecco > 0.05) = nan;
% npp_trends_pct_hycom(npp_trends_sig_hycom > 0.05) = nan;
% par_trends(par_trends_sig > 0.05) = nan;
% ws_trends(ws_trends_sig > 0.05) = nan;
% sla_trends(sla_trends_sig > 0.05) = nan;

for i = 1 : 6
    latIdx = lat < latBand(i, 1) & lat > latBand(i, 2);
    for j = 1 : 4
        %% NEMO
        c1_npp = nan( 120, 3 ); c1_npp_sig = nan( 120, 3 ); % NPP_NEMO
        c1_npp(lat > 0, :) = npp_trends_pct_nemo(lat > 0, mmN(j, :)); c1_npp_sig(lat > 0, :) = npp_trends_sig_nemo(lat > 0, mmN(j, :));
        c1_npp(lat < 0, :) = npp_trends_pct_nemo( lat < 0, mmS(j, :)); c1_npp_sig(lat < 0, :) = npp_trends_sig_nemo( lat < 0, mmS(j, :));
        c1_npp = c1_npp(latIdx, :); c1_npp_sig = c1_npp_sig(latIdx, :); 
        c1_npp = c1_npp(:); c1_npp_sig = c1_npp_sig(:);
        idx_NPP_sig = ~isnan( c1_npp );
        c1_npp = c1_npp( idx_NPP_sig );
        c1_npp_sig = c1_npp_sig( idx_NPP_sig );

        c1_mld = nan( 120, 3 ); c1_mld_sig = nan( 120, 3 ); % MLD_NEMO
        c1_mld(lat > 0, :) = mld_trends_pct_nemo(lat > 0, mmN(j, :)); c1_mld_sig(lat > 0, :) = mld_trends_sig_nemo(lat > 0, mmN(j, :));
        c1_mld(lat < 0, :) = mld_trends_pct_nemo( lat < 0, mmS(j, :)); c1_mld_sig(lat < 0, :) = mld_trends_sig_nemo( lat < 0, mmS(j, :));
        c1_mld = c1_mld(latIdx, :); c1_mld_sig = c1_mld_sig(latIdx, :); 
        c1_mld = c1_mld(:); c1_mld_sig = c1_mld_sig(:);
        c1_mld = c1_mld( idx_NPP_sig );
        c1_mld_sig = c1_mld_sig( idx_NPP_sig );

        %% ECCO
        c2_npp = nan( 120, 3 ); c2_npp_sig = nan( 120, 3 ); % NPP_ECCO
        c2_npp(lat > 0, :) = npp_trends_pct_ecco(lat > 0, mmN(j, :)); c2_npp_sig(lat > 0, :) = npp_trends_sig_ecco(lat > 0, mmN(j, :));
        c2_npp(lat < 0, :) = npp_trends_pct_ecco( lat < 0, mmS(j, :)); c2_npp_sig(lat < 0, :) = npp_trends_sig_ecco( lat < 0, mmS(j, :));
        c2_npp = c2_npp(latIdx, :); c2_npp_sig = c2_npp_sig(latIdx, :); 
        c2_npp = c2_npp(:); c2_npp_sig = c2_npp_sig(:);
        idx_NPP_sig = ~isnan( c2_npp );
        c2_npp = c2_npp( idx_NPP_sig );
        c2_npp_sig = c2_npp_sig( idx_NPP_sig );
        
        c2_mld = nan( 120, 3 ); c2_mld_sig = nan( 120, 3 ); % MLD_ECCO
        c2_mld(lat > 0, :) = mld_trends_pct_ecco(lat > 0, mmN(j, :)); c2_mld_sig(lat > 0, :) = mld_trends_sig_ecco(lat > 0, mmN(j, :));
        c2_mld(lat < 0, :) = mld_trends_pct_ecco( lat < 0, mmS(j, :)); c2_mld_sig(lat < 0, :) = mld_trends_sig_ecco( lat < 0, mmS(j, :));
        c2_mld = c2_mld(latIdx, :); c2_mld_sig = c2_mld_sig(latIdx, :); 
        c2_mld = c2_mld(:); c2_mld_sig = c2_mld_sig(:);
        c2_mld = c2_mld( idx_NPP_sig );
        c2_mld_sig = c2_mld_sig( idx_NPP_sig );

        %% HYCOM
        c3_npp = nan( 120, 3 ); c3_npp_sig = nan( 120, 3 ); % NPP_HYCOM
        c3_npp(lat > 0, :) = npp_trends_pct_hycom(lat > 0, mmN(j, :)); c3_npp_sig(lat > 0, :) = npp_trends_sig_hycom(lat > 0, mmN(j, :));
        c3_npp(lat < 0, :) = npp_trends_pct_hycom( lat < 0, mmS(j, :)); c3_npp_sig(lat < 0, :) = npp_trends_sig_hycom( lat < 0, mmS(j, :));
        c3_npp = c3_npp(latIdx, :); c3_npp_sig = c3_npp_sig(latIdx, :); 
        c3_npp = c3_npp(:); c3_npp_sig = c3_npp_sig(:);
        idx_NPP_sig = (~isnan( c3_npp ));
        c3_npp = c3_npp( idx_NPP_sig );
        c3_npp_sig = c3_npp_sig( idx_NPP_sig );

        c3_mld = nan( 120, 3 ); c3_mld_sig = nan( 120, 3 ); % MLD_HYCOM
        c3_mld(lat > 0, :) = mld_trends_pct_hycom(lat > 0, mmN(j, :)); c3_mld_sig(lat > 0, :) = mld_trends_sig_hycom(lat > 0, mmN(j, :));
        c3_mld(lat < 0, :) = mld_trends_pct_hycom( lat < 0, mmS(j, :)); c3_mld_sig(lat < 0, :) = mld_trends_sig_hycom( lat < 0, mmS(j, :));
        c3_mld = c3_mld(latIdx, :); c3_mld_sig = c3_mld_sig(latIdx, :); 
        c3_mld = c3_mld(:); c3_mld_sig = c3_mld_sig(:);
        c3_mld = c3_mld( idx_NPP_sig );
        c3_mld_sig = c3_mld_sig( idx_NPP_sig );

        %%
        sig_level = 0.05;
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

save( 'senMLD_v2_p005.mat', 'senMLD_NEMO', 'senMLD_ECCO', 'senMLD_HYCOM', 'senPt_nemo', 'senPt_ecco', 'senPt_hycom' )

