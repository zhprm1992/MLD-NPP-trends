% Used for: CbPM NPP (and growth rate, carbon, Chl:C) per 1 deg cell and month, integrated over the euphotic zone, with MLD-NEMO, grid_npp_nemo\npp_YYYY.mat; run once per MLD product (PATH_IN_mld and PATH_OT_npp: NEMO one_d03\mld -> grid_npp_nemo, ECCO MLD_ECCO\one -> grid_npp_ecco, HYCOM MLD_HYCOM\one_gat\mld -> grid_npp_hycom)
% Original: Code\npp_cbpm\grid_npp_1M_1dg_toED_ot_modelMLD.m (data paths updated to I:\Prj; computation unchanged)
% NPP - CbPM monthly using model MLD
% Output:
%   - mu growth rates
%   - NPP
%
% Reference: 
% Westberry et al., (2008), Carbon‐based primary productivity modeling with vertically resolved photoacclimation. Global Biogeochemical Cycles
%
% Haipeng Zhao
% 11/23/2024
% 5/28/2025
% 6/26/2025

clear
clc

PATH_IN_bbp443  = 'I:\Prj\npp_cbpm_1M_1dg\grid_bbp443';
PATH_IN_chl = 'I:\Prj\npp_cbpm_1M_1dg\grid_chl';
PATH_IN_k490 = 'I:\Prj\npp_cbpm_1M_1dg\grid_kd490';
PATH_IN_par = 'I:\Prj\npp_cbpm_1M_1dg\grid_par';
PATH_IN_mld = 'I:\Prj\MLD_1M_NEMO\one_d03\mld'; % MLD PATH 
% PATH_IN_mld = 'I:\Prj\MLD_HYCOM\one_gat\mld';
% PATH_IN_mld = 'I:\Prj\MLD_ECCO\one';

PATH_OT_npp = 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\grid_npp_nemo'; % OUTPUT PATH

parFraction = [0.0029, 0.0032, 0.0035, 0.0037, 0.0037, 0.0036, 0.0032, 0.0030, 0.0024];
r = 0.1;
uMax = 2; % Maximum growth rate (d-1)
LAMBDA = [400, 412, 443, 490, 510, 555, 625, 670, 700];
X = [0.11748, 0.122858, 0.107212, 0.07242, 0.05943, 0.03996, 0.04000, 0.05150, 0.03000];
E = [0.64358, 0.653270, 0.673358, 0.68955, 0.68567, 0.64204, 0.64700, 0.69500, 0.60000];
Kw = [0.01042, 0.007932, 0.009480, 0.01660, 0.03385, 0.06053, 0.28400, 0.43946, 0.62438];

chlCarbonMax = 0; % Max chl:carbon ratio
chlCarbonSat = 0; % Satellite chl:carbon ratio
y0 = 0.0003; % min chl:c when mu = 0
carbon = 0; % bbp converted to carbon
z = 1 : 300;

file_bbp443 = dir( [PATH_IN_bbp443, '\', '*.mat'] ); % bbp443
file_k490 = dir( [PATH_IN_k490, '\', '*.mat'] ); % par
file_par = dir( [PATH_IN_par, '\', '*.mat'] ); % par
file_chl = dir( [PATH_IN_chl, '\', '*.mat'] ); % chl
file_mld = dir( [PATH_IN_mld, '\', '*.mat'] ); % mld

% d = load( 'I:\Prj\WOA23\nitracline\zno3_woa23_05muM.mat' );
% d = load( 'I:\Prj\WOA23\nitracline\zno3_woa23_mod' );
% d = load( 'I:\Prj\WOA23\nitracline\zno3_woa23_Arteaga' );
d = load( 'I:\Prj\WOA23\nitracline\zno3_woa23_3muM.mat' );

zno3 = d.zno3;
lat_zon3 = d.lat;
lon_zon3 = d.lon;
zno3 = [zno3(:, 181 : end, :), zno3(:, 1 : 180, :)];

for i = 1 : 20
    yyyy = num2str( i + 2002 );
    nppG1dg = nan( 180, 360, 'single' ); 
    muG1dg = nan( 180, 360, 'single' ); 
    carbonG1dg = nan( 180, 360, 'single' ); 
    ntFUNG1dg = nan( 180, 360, 'single' ); 
    IgFUNG1dg = nan( 180, 360, 'single' ); 
    chlCarbonG1dg = nan( 180, 360, 'single' ); 
    chlCarbonMxG1dg = nan( 180, 360, 'single' ); 
    chlG1dg = nan( 180, 360, 'single' ); 
    parG1dg =  zeros( 180, 360, 'single' );
    par_mld = zeros( 180, 360, 12 );
    for j = 1 : 12
        mm = num2str( j );

        load( [file_bbp443(i).folder, '\', file_bbp443(i).name], 'bbp443' );
        load( [file_k490(i).folder, '\', file_k490(i).name], 'kd490' );
        load( [file_par(i).folder, '\', file_par(i).name], 'par' );
        load( [file_chl(i).folder, '\', file_chl(i).name], 'chl' );
        load( [file_mld(i).folder, '\', file_mld(i).name], 'mldG1deg' ); % MLD
    
        bbp443 = bbp443(:, :, j);
        kd490 = kd490(:, :, j);
        par = par(:, :, j);
        chl = chl(:, :, j);
        mld = mldG1deg(:, :, j);
        fName = 'I:\Prj\ml_20030101-20030131.nc';
        lat = single( ncread( fName, 'LATITUDE' ) );
        lon = single( ncread( fName, 'LONGITUDE' ) );
        
        Kd = nan( 180, 360, 'single' );
        Kdif = nan( 180, 360, 'single' );
        Klambda = nan( 180, 360, 9, 'single' ); % K(lambda) diffuse attenuation
        Eo = nan( 180, 360, 9, 'single' ); % Spectral PAR
        Ez_mld = nan( 180, 360, 9, 'single' ); %  

        prcnt = zeros( 180, 360, 300, 'single' );
        chl_C = zeros( 180, 360, 300, 'single' );
        mu = zeros( 180, 360, 300, 'single' );
        chl_z = zeros( 180, 360, 300, 'single' );   
        par_z = zeros( 180, 360, 300, 'single' );   
        c_z = zeros( 180, 360, 300, 'single' ); 
        nutTempFun_z =  zeros( 180, 360, 300, 'single' );
        IgFun_z =  zeros( 180, 360, 300, 'single' );
        chlCarbonSat_z =  zeros( 180, 360, 300, 'single' );
        chlCarbonMax_z =  zeros( 180, 360, 300, 'single' );

        Ezlambda = zeros( 180, 360, 9, 300, 'single' ); 
        mzeu = single( ~isnan( par ) ); 
        mzeu(mzeu == 0) = nan;       
        ppz = zeros( 180, 360, 300, 'single' );

        dl = daylength( lat, j, i + 2002 ) * 24;  % in decimal hour

        %% Initialize values
        for k = 1 : length( LAMBDA )
            Klambda(:, :, k) = austinPetzold_1986( LAMBDA(k), kd490 ); % diffuse attenuation coefficient K(lambda)
            Eo(:, :, k) = par .* parFraction(k); % spectral PAR
            Ez_mld(:, :, k) = Eo(:, :, k) .* 0.975 .* exp( -Klambda(:, :, k) .* mld / 2 );
        end          

        for k = 1 : length( LAMBDA ) - 1
            par_mld(:, :, j) = par_mld(:, :, j) + (LAMBDA(k + 1) - LAMBDA(k)) .* (Ez_mld(:, :, k + 1) + Ez_mld(:, :, k)) / 2;
        end
        
        par_mld(:, :, j) = par_mld(:, :, j) ./ dl;
        IgFun = 1 - exp( -5 * par_mld(:, :, j) );

        %% Estimate phytoplankton carbon
        bbp443(bbp443 < 0.00035) = 0.00036; % Original code uses 0.00035
        carbon = 13000 * (bbp443 - 0.00035);
        
        chlCarbonSat = chl ./ carbon;
        chlCarbonSat(chlCarbonSat < y0) = y0;

        chlCarbonMax = 0.022 + (0.045 - 0.022) .* exp( -3.0 .* par_mld(:, :, j) );
        delChlC = chlCarbonMax - chlCarbonSat;
        chlCarbonSat(delChlC < 0) = chlCarbonMax(delChlC < 0); % Original code doesn't have this step

        nutTempFun = (chlCarbonSat - y0) ./ (chlCarbonMax - y0);            
        delChlC(delChlC < 0) = 0; % Original code doesn't have this step
        
        %% Calculate Kd offset carry through to depth non-chl attenuation
        for k = 1 : 9
            Kbio = X(k) * chl .^ E(k);
            Kd(:, :, k) = Kw(k) + Kbio;
            tmp = Klambda(:, :, k) - Kd(:, :, k);
            Kdif(:, :, k) = max( 0, tmp, "includemissing" );            
        end

        %% Integrate down the water column in one-meter steps (mld or 300 m)
        for m = 1 : length( z )
            z_m = z(m);
            aboveMLD = mld;
            aboveMLD(z_m <= mld) = 1;
            aboveMLD(z_m > mld) = 0;

            belowMLD = aboveMLD;
            belowMLD(aboveMLD == 0) = 1;
            belowMLD(aboveMLD == 1) = 0;

            chl_C(:, :, m) = chl_C(:, :, m) + aboveMLD .* chlCarbonSat;
            chl_z(:, :, m) = chl_z(:, :, m) + aboveMLD .* chl_C(:, :, m) .* carbon;
            mu1 = uMax .* nutTempFun .* IgFun;
            mu1 = min( mu1, uMax, 'includenan' );
            mu(:, :, m) = aboveMLD .* mu1;         

            for k = 1 : 9
                Ezlambda(:, :, k, m) = Ezlambda(:, :, k, m) + aboveMLD .* Eo(:, :, k) .* 0.975 .* exp( -Klambda(:, :, k) * z_m );
            end

            c_z(:, :, m) = aboveMLD .* carbon;
            nutTempFun_z(:, :, m) = aboveMLD .* nutTempFun;
            IgFun_z(:, :, m) = aboveMLD .* IgFun;
            chlCarbonSat_z(:, :, m) = aboveMLD .* chlCarbonSat;
            chlCarbonMax_z(:, :, m) = aboveMLD .* chlCarbonMax;

            parz = 0;
            for k = 1 : 8
                parz = parz + (LAMBDA(k + 1) - LAMBDA(k)) * (Ezlambda(:, :, k + 1, m) + Ezlambda(:, :, k, m)) / 2;
            end
            par_z(:, :, m) = par_z(:, :, m) + aboveMLD .* parz;

            if m > 1
                chlz_prev = chl_z(:, :, m - 1);
                for k = 1 : 9
                    Kbio = X(k) .* chlz_prev .^ E(k);
                    Kd(:, :, k) = Kw(k) + Kbio;
                    Kd(:, :, k) = Kd(:, :, k) + Kdif(:, :, k);    % NOTE: Here Kdif (non-chl kd) is assumted constant vertically
                    Ezlambda(:, :, k, m) = Ezlambda(:, :, k, m) + belowMLD .* Ezlambda(:, :, k, m - 1) .* exp( -Kd(:, :, k) * 1 );
                end
                
                parz = 0;
                for k = 1 : 8
                    parz = parz + (LAMBDA(k + 1) - LAMBDA(k)) * (Ezlambda(:, :, k + 1, m) + Ezlambda(:, :, k, m)) / 2;
                end
                par_z(:, :, m) = par_z(:, :, m) + belowMLD .* parz;

                deltaZ = max( 0, zno3(:, :, j) - z_m );
                chl_C(:, :, m) = chl_C(:, :, m) + belowMLD .* (0.022 + (0.045 - 0.022) .* exp( -3.0 .* par_z(:, :, m) ./ dl )  );
                ch_C_m = chl_C(:, :, m) - belowMLD .* (delChlC .* (1 - exp( -0.075 * deltaZ )));
                ch_C_m(ch_C_m < y0) = y0;
                chl_C(:, :, m) = ch_C_m;

                IgFun_z(:, :, m) = IgFun_z(:, :, m) +  belowMLD .* (1 - exp( -5.0 * par_z(:, :, m) ./ dl ));
                chlCarbonMax_z(:, :, m) = chlCarbonMax_z(:, :, m) + belowMLD .* (0.022 + (0.045 - 0.022) .* exp( -3.0 .* par_z(:, :, m) ./ dl ));
                
                nutTempFun_belowMLD = ((chl_C(:, :, m) - y0) ./ (chlCarbonMax_z(:, :, m) - y0));
                nutTempFun_z(:, :, m) = nutTempFun_z(:, :, m) + belowMLD .* nutTempFun_belowMLD; % NOTE: Orignal code does include this step.

                % mask = (mld + 5) > zno3(:, :, j); 
                % tmp = nutTempFun_z(:, :, m);     
                % tmp(mask) = nutTempFun(mask);     
                % nutTempFun_z(:, :, m) = tmp;  

                mu_ = uMax .* nutTempFun_z(:, :, m) .* IgFun_z(:, :, m);  
                mu_ = min( mu_, uMax, 'includenan' );       
                mu(:, :, m) = mu(:, :, m) + belowMLD .* mu_;
        
                mu_prev = mu(:, :, m - 1);
                c_z_below = carbon .* (mu_prev ./ r);
                c_z_below(mu_prev >= r) = carbon(mu_prev >= r);
                c_z(:, :, m) = c_z(:, :, m) + belowMLD .* c_z_below;
        
                chl_z(:,:,m) = chl_z(:,:,m) + belowMLD .* chl_C(:,:,m) .* c_z(:,:, m);
            end

            prcnt(:, :, m) = par_z(:, :, m) ./ (par .* 0.975);

            % mzeu update (depth where prcnt drops below 0.01)
            is_new = prcnt(:, :, m) >= 0.01;
            mzeu(is_new) = m;
        
            % z_eu interpolation
%             if m > 1
%                 interp_mask = (z_eu == 0) & (prcnt(:, :, m) < 0.01) & (mzeu > 0);
%                 if any(interp_mask, 'all')
%                     p0 = prcnt(:, :, mzeu);
%                     p1 = prcnt(:, :, m);
%                     z0 = z(mzeu);
%                     z1 = z(m);
%                     frac = (p0 - 0.01) ./ (p0 - p1);
%                     z_eu(interp_mask) = z0 + (z1 - z0) .* frac;
%                 end
%             end
        
            % ppz
            ppz(:, :, m) = mu(:, :, m) .* c_z(:, :, m);
        end

        % c_z( c_z == 0 ) = nan;
        % mu( c_z == 0 ) = nan;
        % nutTempFun_z( c_z == 0 ) = nan;
        % IgFun_z( c_z == 0 ) = nan;
        % chlCarbonSat_z( c_z == 0 ) = nan;
        % chlCarbonMax_z( c_z == 0 ) = nan;

        %% Calculate NP
        mask = mzeu < 186;
        nppG1dg(:, :, j) = mask .* trapz( 1 : 300, ppz, 3 );
        carbonG1dg(:, :, j) = mask .* trapz( 1 : 300, c_z, 3 );
        muG1dg(:, :, j) = mask .* nanmean( mu, 3 );
        ntFUNG1dg(:, :, j) = mask .* nanmean( nutTempFun_z, 3 );
        IgFUNG1dg(:, :, j) = mask .* nanmean( IgFun_z, 3 );
        chlCarbonG1dg(:, :, j) = mask .* nanmean( chlCarbonSat_z, 3 );
        chlCarbonMxG1dg(:, :, j) = mask .* nanmean( chlCarbonMax_z, 3 );
        chlG1dg(:, :, j) = mask .* nanmean( chl_z, 3 );
        parG1dg(:, :, j) = mask .* nanmean( par_z, 3 );

        nppG1dg( nppG1dg == 0 ) = nan;
        muG1dg( muG1dg == 0 ) = nan;
    end
    nppG1dg( nppG1dg == 0 ) = nan;
    par_mld( par_mld == 0) = nan;
    filename_ot = ['npp_', yyyy]; 
    save( fullfile( [PATH_OT_npp, '\', filename_ot] ), 'muG1dg', 'carbonG1dg', 'chlCarbonMxG1dg', 'nppG1dg',...
        'ntFUNG1dg', 'IgFUNG1dg', 'chlCarbonG1dg', 'par_mld', 'chlG1dg', 'parG1dg', 'lon', 'lat', '-v7.3'  ) % Grid = [180, 360], Lat: (-80 , 90)
end

%%
function DL = daylength( lat, mm, yyyy )
% lat: Latitude
% mm: Month
% yyyy: Year
    PI = 3.141592653589793;

    num_days = eomday( yyyy, mm );
    dl = nan( 180, 360, num_days );

    nDays = yeardays( yyyy ); % e.g. 365 
    lat = repmat( lat, 1, 360 );
    for i = 1 : num_days
        yd = day( datetime( yyyy, mm, i ), 'dayofyear' );

        % Declination (radians)
        decl = 23.5 * cos(2 * PI * (yd - 172) / nDays) * PI / 180;

        % Daylength (not used below, but preserved if needed)
        DL = -1 * tan( lat * PI / 180 ) * tan( decl );
        DL = min( max( DL, -1 ), 1 );
        dl(:, :, i) = acos( DL ) / PI; % Fraction of day
    end

    DL = nanmean( dl, 3 );
end