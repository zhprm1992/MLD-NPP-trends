% Used for: CAFE NPP per 1 deg cell and month, integrated over the euphotic zone; set here for MLD-ECCO, run once per MLD product (PATH_IN_mld, PATH_OT_npp). Line 266 is the version that produced the stored CAFE grids (verified: June 2010 NEMO and ECCO reproduced exactly); the later 11/4/2026 change of that line (E_wz -> E_wz2) was not used for any published result
% Original: Code\npp_cafe\grid_npp_1M_1dg_toED_ot_modelMLD.m (data paths updated to I:\Prj; computation unchanged)
% Output CAFE derived NPP. This version uses model MLD 
% Haipeng Zhao
% 1/6/2024
% 5/21/2025   NPP to euphotic zone
% 5/2/2026 Fix some bugs
% 11/4/2026 line264: E_wz -> E_wz2
%
% Parameter related to eupohotic zone
% decl = 0;     % Solar declination [rads]
% m0 = 0;       % Coefficient to calculate kd (Lee et al. 2005) [m-1]
% DL = 0;        % Daylength [day]
% solzen = 0;          % Solar zenith angle [rads]
% kd = zeros(31, 1); % Downwelling attenuation coefficient of irradiance [m-1]
% zeu = 0;                    % Euphotic depth [m]
% kdpar = 0;                % Downwelling attenuation coefficient of PAR [m-1]
%
% Paramters related to light and absorbed photons through depth, time, and wavelength
% tseq = zeros( 101, 1 );         % Fractional time of day
% zseq = zeros( 101, 1 );         % Euphotic depth divided into 101 increments
% delz = 0;                                    % Depth of zseq (i.e. zeu/101)
% absorbed_photons = 0;         % Absorbed photons calculated  in Step 3
% PAR_noon = zeros( 31, 1 );  % PAR at solar noon[wavelength] [mol photons m-2 day-1 wv-]
% AP_tzw = zeros( 101, 101, 31 );     % Absorbed photons[depth, time, wavelength]
% AP_tz = zeros( 101, 101 );                % Absorbed photons [depth, time]
% AP_tz2 = zeros( 101, 101 );                    % Absorbed photons[depth, time]
% AP_z = zeros( 101, 1 );           % Absorbed photons[depth]
% AP_z2 = zeros(101, 1);           % Absorbed photons[depth]
% AP = 0;                                         % Absorbed photons [mol m-2 day-1]
% E_tzw = zeros( 101, 101, 31 );     % Irradiance [depth, time, wavelength]
% E_tz = zeros( 101, 101 );                % Irradiance [depth, time]
% Eu = 0;                                                  % conversion factor (Eu) that converts downwelling planar irradiance to scalar irradiance
%
% Parameters related to derivation of Ek
% Ek = zeros( 101, 1 );  % Photoacclimation parameter[depth] (Behrenfeld et al. 2016) [mol m-2 day-1]
% IML = 0;        % Median irradiance in the mixed layer [umol m-2 day-1]
% Eg = 0;           % Growth irradiance at a specific depth
% Eg_mld;          % Growth irradiance at a specific depth
%
% Parameters used to convert Ek into the spectrally explicit equivalent KPUR
% KPUR = zeros( 101, 1 );        % Spectrally explicit Ek
% numerator;                                % Numerator to calculate spectral correction factor
% denominator;                            % Denominator to calculate spectral correction factor
% mean_aphi = 0;                         % Spectrally averaged absorption due to phytoplankton
%
% Parameters used to calculate the maximum quantum yield of net carbon fixation
% phimax = zeros( 101, 1 );                                   % maximum quantum yield of net carbon fixation [mol C (mol photons)-1]
% phirange = [0.018, 0.030]';                                 % range of phimax values  
% Ekrange = [150 * 0.086400, 10 * 0.0864]';    % relates phimax to Ek
% slope = 0;                                                                  % relates phimax to Ek
% aphi_fact = zeros( 101, 1 );    %Parameter to increase aphi beneath the MLD
%
% Net primary production [mg C m-2 day-1] */
% NPP_tz = zeros( 101, 101 );
% NPP_z = zeros(101, 1 );
% NPP = 0;

clear
clc

%% Step 1: Declare local variables
PI = 3.141592653589793;
num_z = 101;
num_t = 51;
num_w = 31;

% Wavelenght
wv = [400, 410, 420, 430, 440, 450, 460, 470, 480, 490, 500, 510, 520, 530, 540, 550, 560, ...
		 570, 580, 590, 600, 610, 620, 630, 640, 650, 660, 670, 680, 690, 700]';

% Absorptioin due to pure water
aw = [0.00663, 0.00473, 0.00454, 0.00495, 0.00635, 0.00922, 0.00979, 0.0106, 0.0127,... 
		  0.015, 0.0204, 0.0325, 0.0409, 0.0434, 0.0474, 0.0565, 0.0619, 0.0695,... 
		  0.0896, 0.1351, 0.2224, 0.2644, 0.2755, 0.2916, 0.3108, 0.34, 0.41, 0.439, 0.465,... 
		  0.516, 0.624]';

% Spectral shape of aphi.  aphi = A * chl ^ E
A_Bricaud = [0.0241, 0.0287, 0.0328, 0.0359, 0.0378, 0.0350, 0.0328, 0.0309, 0.0281,...
			 0.0254, 0.0210, 0.0162, 0.0126, 0.0103, 0.0085, 0.0070, 0.0057, 0.0050, ...
			 0.0051, 0.0054, 0.0052, 0.0055, 0.0061, 0.0066, 0.0071, 0.0078, 0.0108, ...
			 0.0174, 0.0161, 0.0069, 0.0025]';

E_Bricaud = [0.6877, 0.6834, 0.6664, 0.6478, 0.6266, 0.5993, 0.5961, 0.5970, 0.5890,... 
			 0.6074, 0.6529, 0.7212, 0.7939, 0.8500, 0.9036, 0.9312, 0.9345, 0.9298, ...
			 0.8933, 0.8589, 0.8410, 0.8548, 0.8704, 0.8638, 0.8524, 0.8155, 0.8233, ...
			 0.8138, 0.8284, 0.9255, 1.0286]';

% Spectral shape of PAR
 PAR_spectrum = [0.00227, 0.00218, 0.00239, 0.00189, 0.00297, 0.00348, 0.00345, ...
			    0.00344, 0.00373, 0.00377, 0.00362, 0.00364, 0.00360, 0.00367, ...
			    0.00354, 0.00368, 0.00354, 0.00357, 0.00363, 0.00332, 0.00358, ...
			    0.00357, 0.00359, 0.00340, 0.00350, 0.00332, 0.00342, 0.00347, ...
			    0.00342, 0.00290, 0.00314]';

phirange = [0.018, 0.030]';                                 % range of phimax values  
Ekrange = [150 * 0.086400, 10 * 0.0864]';    % relates phimax to Ek

PATH_IN_sst = 'I:\Prj\npp_cafe_1M_1dg\grid_sst_1M_1dg_mat';    % SET UP PATH OF INPUT
PATH_IN_aph443 = 'I:\Prj\npp_cafe_1M_1dg\grid_aph443_1M_1dg_mat';
PATH_IN_adg443 = 'I:\Prj\npp_cafe_1M_1dg\grid_adg443_1M_1dg_mat';
PATH_IN_bbp443 = 'I:\Prj\npp_cafe_1M_1dg\grid_bbp443_1M_1dg_mat';
PATH_IN_PAR = 'I:\Prj\npp_cafe_1M_1dg\grid_par_1M_1dg_mat';
PATH_IN_chl = 'I:\Prj\npp_cafe_1M_1dg\grid_chl_1M_1dg_mat';
PATH_IN_bbps = 'I:\Prj\npp_cafe_1M_1dg\grid_bbps_1M_1dg_mat';
% PATH_IN_mld = 'I:\Prj\mld_1M_NEMO\one_d03\mld';
% PATH_IN_mld = 'I:\Prj\MLD_HYCOM\one_gat\mld';
PATH_IN_mld = 'I:\Prj\MLD_ECCO\one';

PATH_OT_npp = 'I:\Prj\npp_cafe_1M_1dg\grid_npp_ECCO_1M_1dg_mat'; % SET UP PATH OF OUTPUT

file_aph443 = dir( [PATH_IN_aph443, '\', '*.mat'] ); % aph443
file_adg443 = dir( [PATH_IN_adg443, '\', '*.mat'] ); % adg443
file_bbp443 = dir( [PATH_IN_bbp443, '\', '*.mat'] ); % bbp443
file_bbps = dir( [PATH_IN_bbps, '\', '*.mat'] ); % bbps
file_par = dir( [PATH_IN_PAR, '\', '*.mat'] ); % par
file_chl = dir( [PATH_IN_chl, '\', '*.mat'] ); % chl
file_sst = dir( [PATH_IN_sst, '\', '*.mat'] ); % sst
file_mld = dir( [PATH_IN_mld, '\', '*.mat'] ); % mld

for i = 1 : 20 % 2003 - 2022
    yyyy = num2str( i + 2002 );
    for j = 1 : 12 
        mm = num2str( j );
        % Load data
        load( [file_aph443(i).folder, '\', file_aph443(i).name], 'aph443' );
        load( [file_adg443(i).folder, '\', file_adg443(i).name], 'adg443' );
        load( [file_bbp443(i).folder, '\', file_bbp443(i).name], 'bbp443' );
        load( [file_bbps(i).folder, '\', file_bbps(i).name], 'bbps' );
        load( [file_par(i).folder, '\', file_par(i).name], 'par' );
        load( [file_chl(i).folder, '\', file_chl(i).name], 'chl' );
        load( [file_sst(i).folder, '\', file_sst(i).name], 'sst' );
        load( [file_mld(i).folder, '\', file_mld(i).name], 'mldG1deg' ); % MLD    
    
        aph443 = aph443(:, :, j);
        adg443 = adg443(:, :, j);
        bbp443 = bbp443(:, :, j);
        bbps = bbps(:, :, j);
        sst = sst(:, :, j);
        par = par(:, :, j);
        chl = chl(:, :, j);
        mld = mldG1deg(:, :, j);
        fName = 'I:\Prj\MLD_ARGO\one\0.03sigmaTheta_ml_2003.nc\ml_20030101-20030131.nc';
        lat = single( ncread( fName, 'LATITUDE' ) );
        
        [r, c] = size(chl);
        bbw = nan( r, c, 31, 'single' );
        aphi = nan( r, c, 31, 'single' );
        a = nan( r, c, 31, 'single' );
        bb = nan( r, c, 31, 'single' );
        
        %% Step: 2 & 3
        for w = 1 : numel( wv )
            % 2. Derive IOPs at 10 nm increments from 400 to 700 nm. Comment: Currently assume salinity is 32.5 ppmil
            bbw(:, :, w) = betasw_ZHH2009( wv(w), 32.5, sst ) ./ 2;
            aphi(:, :, w) = aph443 .* A_Bricaud(w) .* (chl .^ E_Bricaud(w)) ./ (0.03711 * chl .^ 0.61479);
            a(:, :, w) = aw(w) + aphi(:, :, w) + adg443 .* exp( -0.018 * (wv(w) - 443) );
            bb(:, :, w) = bbw(:, :, w) + bbp443 .* (443 / wv(w)) .^ bbps;
        end
    
        % 3. Calculate Absorbed Energy (CAVEAT: THE CODE INTEGRATES FRACTION OF PHOTONS which may need a correction)
        aphi_over_a = aphi ./ a .* reshape( PAR_spectrum, 1, 1, [] );
        absorbed_photons = 10 * trapz( 1 : 31, aphi_over_a, 3 );
        absorbed_photons = absorbed_photons .* par * 0.95;
    
        %% 4. Derive Kd following Lee et al 2005, Eq. 11 and kdpar from Morel et al 2007, Eq. 9
        [kd, kdpar, zeu, DL, m0] = computeMonthlyOptics( lat, a, bb, par, j, i + 2002 );
    
        %% 5. Derive conversion factor (Eu) that converts downwelling planar irradiance to scalar irradiance 
        zeu_ceil = ceil(zeu);               
        z_idx = reshape( 0 : (num_z - 1), 1, 1, [] ); 
        zseq = zeu_ceil ./ (num_z - 1) .* z_idx;
        delz = zseq(:, :, 2) - zseq(:, :, 1);          
                 
        E_wz = nan( r, c, 31, num_z, 'single' );
        AP_wz = nan( r, c, 31, num_z, 'single' );
    
        for z = 1 : size( zseq, 3 )
            z_slice = zseq(:, :, z);
            z_slice = repmat( z_slice, 1, 1, 31 );
            par_noon = PI / 2 * par * 0.95 .* reshape( PAR_spectrum, 1, 1, [] );
            E_wz(:, :, :, z) = 2 / PI * par_noon .* exp( -1 * kd .* z_slice ); % [lat, lon, wavelength, depth]
            AP_wz(:, :, :, z) = E_wz(:, :, :, z) .* aphi;
        end
    
        % integrate AP through wavelength
        AP_z = squeeze( 10 * trapz( 1 : 31, AP_wz, 3 ) );
        E_z = squeeze( 10 * trapz( 1 : 31, E_wz, 3 ) );
    
        % INtegrate AP through depth
        AP = delz .* trapz( 1 : 101, AP_z, 3 );
    
        % Derive upwelling irradiance
        Eu = absorbed_photons ./ AP;
        Eu = repmat( Eu, 1, 1, 101 );
        % Multiply E_z by Eu
        
        E_z = E_z .* Eu;
        AP_z = AP_z .* Eu;
    
        %% 6. CALCULATE EK through depth (modified from Behrenfeld et al. 2016)
        IML = (par * 0.95 ./ (DL * 24)) .* exp( -0.5 * kdpar .* mld ); % Median irradiance in the mixed layer [umol m-2 day-1]
        
        Ek = 19 .* exp( 0.038 .* ((par .* 0.95 ./ (DL .* 24)) .^ 0.45) ./ kdpar ); 
        Ek = max( Ek, 10 ); % Enforce Ek >= 10 at all depths
        Ek = repmat( Ek, 1, 1, num_z );
    
        iml_term = (1 + exp(-0.15 .* (par .* 0.95 ./ (DL .* 24)))) ./ (1 + exp(-3 .* IML)); 
        iml_term = repmat(iml_term, 1, 1, 101);  
        Ek_update = Ek .* iml_term;  
    
        mask_mld_zeu = mld < zeu;  
        mask_mld_zeu = repmat(mask_mld_zeu, 1, 1, 101);
        Ek(mask_mld_zeu) = Ek_update(mask_mld_zeu);
    
        Eg_mld = (par ./ DL) .* exp( -1 * kdpar .* mld );
        Eg_mld = repmat( Eg_mld, 1, 1, 101 );  
        Eg = (par ./ DL) .* exp( -1 * kdpar .* zseq ); 
        Ek_update = 10 + (Ek - 10) ./ (Eg_mld - 0.1) .* (Eg - 0.1);
    
        z_gt_mld = zseq > repmat(mld, 1, 1, 101);       
        Ek(mask_mld_zeu & z_gt_mld) = Ek_update(mask_mld_zeu & z_gt_mld);
    
        Ek = max( Ek, 10 );
        Ek = Ek * 0.0864;
    
        %% 7. Make Ek spectrally explicit (KPUR) using a spectral correction factor
        mean_aphi = mean( aphi, 3, 'omitnan' );
        E_noon_wz = nan( r, c, 31, num_z, 'single' );
        AP_noon_wz = nan( r, c, 31, num_z, 'single' );
        for z = 1 : num_z
            E_noon_wz(:, :, :, z) = par_noon .* 1 .* exp( -1 * kd .* zseq(:, :, z) );
            AP_noon_wz(:, :, :, z) = E_noon_wz(:, :, :, z) .* aphi;
        end
        
        denominator = squeeze( 10 * trapz( 1 : 31, E_noon_wz, 3 ) ); % Denominator for each depth
        numerator = squeeze( 10 * trapz( 1 : 31, AP_noon_wz, 3 ) ); % Numerator for each depth
       
        kpur = Ek ./ (numerator ./ (denominator .* mean_aphi) / 1.3);
    
        %% 8. Calculate phimax as a function of Ek
        slope = (phirange(2) - phirange(1)) / (Ekrange(2) - Ekrange(1));
        
        phimax = phirange(2) + (Ek - Ekrange(2)) * slope;
        
        phimax(phimax > phirange(2)) = phirange(2);
        phimax(phimax < phirange(1)) = phirange(1); 
            
        %% 9. Scale aphi beneath the MLD to Ek
        aphi_fact = ones( r, c, num_z, 'single' );
        aphi_fact_update = 1 + Ek(:, :, 1) ./ Ek * 0.15;
        mask = mask_mld_zeu & z_gt_mld;
        aphi_fact(mask) =  aphi_fact_update(mask);
    
        wv_3d = reshape( wv, 1, 1, 31 );
        aw_3d = reshape( aw, 1, 1, 31 );
    
        E_wz2 = nan( r, c, num_w, num_z, 'single' );
        AP_wz2 = nan( r, c, num_w, num_z, 'single' );
    
        E_wz2(:, :, :, 1) = E_wz(:, :, :, 1);
        AP_wz2(:, :, :, 1) = E_wz(:, :, :, 1) .* aphi;
        for z = 2 : num_z
            a2 = aw_3d + aphi .* aphi_fact(:, :, z) + adg443 .* exp( -0.018 * (wv_3d - 443) );
            kd2 = m0 .* a2 + 4.18 * (1 - 0.52 * exp( -10.8 * a2 )) .* bb;
            E_wz2(:, :, :, z) = E_wz(:, :, :, z - 1) .* exp( -1 * kd2 .* delz );   % as run for the published CAFE NPP (the 11/4/2026 change E_wz -> E_wz2 was made after the runs)
            AP_wz2(:, :, :, z) = E_wz2(:, :, :, z)  .* aphi .* aphi_fact(:, :, z);
        end
        
        AP_z2 = squeeze( 10 * trapz( 1 : 31, AP_wz2, 3 ) );
        AP_z2 = AP_z2 .* Eu;
    
        E_z2 = squeeze( 10 * trapz( 1 : 31, E_wz2, 3 ) ); % This is not applied in original code (http://orca.science.oregonstate.edu/cafe.code.php)
        E_z2 = E_z2 .* Eu;
    
%         AP_z(mask_mld_zeu) = AP_z_update(mask_mld_zeu);
%         E_z(mask_mld_zeu) = E_z_update(mask_mld_zeu);  
%     
        %% 10. Calculate NPP
        npp_z = phimax .* AP_z2 .* tanh( kpur ./ E_z2 ) * 12000;
        npp =  delz .* trapz( 1 : 101, npp_z, 3 );
        nppG1dg(:, :, j) = npp;
        
        ap = delz .* trapz( 1 : 101, AP_z2, 3 );
        apG1dg(:, :, j) = ap;

    end
    %% Output npp
    filename_ot = ['npp_', yyyy]; 
    save( fullfile( PATH_OT_npp, filename_ot ), 'nppG1dg', 'apG1dg', '-v7.3' )

end


%%
function [bsw, beta90sw]= betasw_ZHH2009( lambda, S, Tc )
% Scattering by pure seawater: Effect of salinity Xiaodong Zhang, Lianbo Hu, 
% and Ming-Xia He, Optics Express, 2009
% lambda (nm): wavelength
% Tc: temperature in degree Celsius, must be a scalar
% S: salinity, must be scalar
% delta: depolarization ratio, if not provided, default = 0.039 will be used.
% betasw: volume scattering at angles defined by theta. Its size is [x y],
% where x is the number of angles (x = length(theta)) and y is the number
% of wavelengths in lambda (y = length(lambda))
% beta90sw: volume scattering at 90 degree. Its size is [1 y]
% bw: total scattering coefficient. Its size is [1 y]
% for backscattering coefficients, divide total scattering by 2

Na = 6.02214076e23;        % Avogadro's constant
Kbz = 1.380649e-23;         % Boltzmann constant
Tk = Tc + 273.15;                % Absolute tempearture
delta = 0.039;
M0 = 18e-3;                           % Molecular weigtht of water in kg/mol

% Refractive index of seawater.  Quan and Fry (1994, Applied Optics)
n_air = 1.0 + (5792105.0 ./ (238.0185 - 1 ./ (lambda / 1e3) .^2) + ...
               167917.0 ./ (57.362 - 1 ./ (lambda / 1e3) .^ 2)) / 1e8;

% Refractive index of seawater. Quan an Fry (1994, Applied Optics)
n = [1.31405, 1.779e-4, -1.05e-6, 1.6e-8, -2.02e-6, 15.868, 0.01155, -0.00423,  -4382, 1.1455e6]';

% Pure seawater
nsw = n(1) + (n(2) + n(3) * Tc + n(4) * Tc.^2) .* S + n(5) * Tc.^2 + ...
      (n(6) + n(7) * S + n(8) * Tc) ./ lambda + n(9) ./ lambda.^2 + n(10) ./ lambda.^3;
nsw = nsw .* n_air;
dnswds = (n(2) + n(3) * Tc + n(4) * Tc.^2 + n(7) ./ lambda) .* n_air;

% Pure water secant bulk Millero (1980, Deep-sea Research) */
kw = 19652.21 + 148.4206 * Tc - 2.327105 * Tc.^2 + 1.360477e-2 * Tc.^3 - ...
     5.155288e-5 * Tc.^4;
Btw_cal = 1 ./ kw;

% Isothermal compressibility from Kell sound measurement in pure water
Btw = (50.88630 + 0.717582 * Tc + 0.7819867e-3 * Tc .^ 2 + 31.62214e-6 * Tc .^ 3 - ...
       0.1323594e-6 * Tc .^ 4 + 0.634575e-9 * Tc .^ 5) ./ (1 + 21.65928e-3 * Tc) * 1e-6;

% Seawater secant bulk
g0 = 54.6746 - 0.603459 * Tc + 1.09987e-2 * Tc .^ 2 - 6.167e-5 * Tc .^ 3;
g1 = 7.944e-2 + 1.6483e-2 * Tc - 5.3009e-4 * Tc .^ 2;
Ks = kw + g0 .* S + g1 .* S .^ 1.5;
  
% Calculate seawater isothermal compressibility from the secant bulk */
IsoComp = (1 ./ Ks) * 1e-5;  % unit is pa

a = [8.24493e-1, -4.0899e-3, 7.6438e-5, -8.2467e-7, 5.3875e-9, ...
     -5.72466e-3, 1.02270e-4, -1.6546e-6, 4.8314e-4];

b = [999.842594, 6.793952e-2, -9.09529e-3, 1.001685e-4, ...
     -1.120083e-6, 6.536332e-9];

% Density for pure water
density = b(1) + b(2) * Tc + b(3) * Tc .^ 2 + b(4) * Tc .^ 3 + b(5) * Tc .^ 4 + b(6) * Tc .^ 5;

% Density for pure seawater
density = density + ((a(1) + a(2) * Tc + a(3) * Tc .^ 2 + a(4) * Tc.^3 + a(5) * Tc .^ 4) .* S + ...
                     (a(6) + a(7) * Tc + a(8) * Tc .^ 2) .* S .^ 1.5 + ...
                     a(9) * S .^ 2);

% dlnawds is partial derivative of natural logarithm of water activity  w.r.t.salinity
dlnawds = (-5.58651e-4 + 2.40452e-7 * Tc - 3.12165e-9 * Tc .^ 2 + 2.40808e-11 * Tc .^ 3) + ...
          1.5 * (1.79613e-5 - 9.9422e-8 * Tc + 2.08919e-9 * Tc .^ 2 - 1.39872e-11 * Tc .^ 3) .* S .^ 0.5 + ...
          2 * (-2.31065e-6 - 1.37674e-9 * Tc - 1.93316e-11 * Tc .^ 2) .* S;

% Density derivative of refractive index from PMH model
n_wat2 = nsw .^ 2; % Square of nsw
base = (nsw / 3) - (1 / 3) ./ nsw; % Base term calculation
DFRI = (n_wat2 - 1) .* (1 + (2 / 3) * (n_wat2 + 2) .* base.^2);

% Volume scattering at 90 degree due to the density fluctuation
beta_df = pi ^ 2 / 2 * (lambda * 1e-9) .^ (-4) * Kbz * Tk .* IsoComp .* DFRI .^ 2 * (6 + 6 * delta) ./ (6 - 7 * delta);

% volume scattering at 90 degree due to the concentration fluctuation
flu_con = S * M0 * (dnswds .^ 2) ./ density ./ (-1 * dlnawds) / Na;
beta_cf = 2 * pi^2 * (lambda * 1e-9) .^ (-4) .* nsw .^ 2 .* flu_con * (6 + 6 * delta) / (6 - 7 * delta);

% Total volume scattering at 90 degree
beta90sw = beta_df + beta_cf;
bsw = 8 * pi / 3 * beta90sw * (2 + delta) / (1 + delta);

end


%%
function [kd, kdpar, zeu, DL, m0] = computeMonthlyOptics( lat, a, bb, PAR, mm, yyyy )
% lat: Latitude
% a: Total absorption coefficient
% bb: Total backscattering
% PAR: Light radiantion
% mm: Month
% yyyy: Year
    PI = 3.141592653589793;

    num_days = eomday( yyyy, mm );
    m0 = nan( 180, 360, num_days );
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

        % Solar zenith angle at sunrise/sunset (degrees)
        solzen = 90 - asin( sin( lat * PI / 180 ) * sin( decl ) - ...
                           cos( lat * PI / 180 ) * cos( decl ) * cos( PI ) ) * 180 / PI;

        % Coefficient to calculate kd [m-1]
        m0(:, :, i) = sqrt( (1 + 0.005 * solzen)  .^ 2 );
    end

    DL = nanmean( dl, 3 );
    m0 = nanmean( m0, 3 ); % monthly average
    
    m0 = repmat( m0, 1, 1, size( a, 3 ) );
    kd = m0 .* a + 4.18 * (1 - 0.52 * exp( -10.8 * a ))  .* bb; % Diffuse attenuation

    kd10 = kd(:, :, 10);
    kdpar = 0.0665 + 0.874 * kd10 - (0.00121 ./ kd10); 

    % Euphotic depth
    zeu = -1 * log(0.1 ./ (PAR * 0.95)) ./ kdpar;

end











