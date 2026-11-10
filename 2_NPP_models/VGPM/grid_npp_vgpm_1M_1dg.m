% Used for: VGPM NPP per 1 deg cell and month (no MLD input), grid_npp_vgpm_1M_1dg_mat
% Original: Code\npp_vgpm\grid_npp_vgpm_1M_1dg.m (data paths updated to I:\Prj; computation unchanged)
% Output standard VGPM derived NPP.
% Haipeng Zhao
% 7/11/2026
%
% Converted from vgpm_c_version.txt (opp_befa), matching
% https://orca.science.oregonstate.edu/vgpm.code.php:
% Pb_opt uses the Behrenfeld & Falkowski (1997) clamped polynomial
% temperature fit. z_eu uses the Morel 2007 model (O'Malley 02.13.2023
% revision), as saved in vgpm_c_version.txt.
%
% opp_befa( chl, irr, sst, dayL ):
%   chl  - Chlorophyll_a surface concentration [mg chl m-3]
%   irr  - Photosynthetically available radiation [Einstein m-2 day-1]
%   sst  - Sea surface temperature [deg C]
%   dayL - Day length [decimal hours]
%   npp  - Net primary production [mg C m-2 day-1]

clear
clc

%% Step 1: Set up paths
PATH_IN_sst = 'I:\Prj\npp_cafe_1M_1dg\grid_sst_1M_1dg_mat';    % SET UP PATH OF INPUT
PATH_IN_PAR = 'I:\Prj\npp_cafe_1M_1dg\grid_par_1M_1dg_mat';
PATH_IN_chl = 'I:\Prj\npp_cafe_1M_1dg\grid_chl_1M_1dg_mat';

PATH_OT_npp = 'I:\Prj\npp_vgpm_1M_1dg\grid_npp_vgpm_1M_1dg_mat'; % SET UP PATH OF OUTPUT

file_par = dir( [PATH_IN_PAR, '\', '*.mat'] ); % par
file_chl = dir( [PATH_IN_chl, '\', '*.mat'] ); % chl
file_sst = dir( [PATH_IN_sst, '\', '*.mat'] ); % sst

for i = 1 : numel( file_chl ) % 2003 - 2022
    yyyy = num2str( i + 2002 );

    % Load data
    load( [file_par(i).folder, '\', file_par(i).name], 'par' );
    load( [file_chl(i).folder, '\', file_chl(i).name], 'chl', 'lat', 'lon' );
    load( [file_sst(i).folder, '\', file_sst(i).name], 'sst' );

    [r, c] = size( chl, 1 : 2 );
    nppG1dg = nan( r, c, 12, 'single' );

    for j = 1 : 12
        chl_m = chl(:, :, j);
        par_m = par(:, :, j);
        sst_m = sst(:, :, j);

        % Day length [decimal hours], monthly mean
        dayL = computeDayLength( lat, r, c, j, i + 2002 );

        %% opp_befa
        % Euphotic depth (Morel 2007 model)
        X = log10( chl_m );
        log10z_eu = 1.524 - 0.436 * X - 0.0145 * X .^ 2 + 0.0186 * X .^ 3;
        z_eu = 10 .^ log10z_eu;

        % Pb_opt from SST (Behrenfeld & Falkowski 1997 polynomial, clamped)
        pb_opt = 1.2956 + 2.749e-1 * sst_m + 6.17e-2 * sst_m .^ 2 - 2.05e-2 * sst_m .^ 3 ...
            + 2.462e-3 * sst_m .^ 4 - 1.348e-4 * sst_m .^ 5 + 3.4132e-6 * sst_m .^ 6 ...
            - 3.27e-8 * sst_m .^ 7;
        pb_opt( sst_m < -10 ) = 0;
        pb_opt( sst_m >= -10 & sst_m < -1 ) = 1.13;
        pb_opt( sst_m > 28.5 ) = 4.00;

        % Irradiance function
        irrFunc = 0.66125 * par_m ./ (par_m + 4.1);

        % Net primary production [mg C m-2 day-1]
        npp = pb_opt .* chl_m .* dayL .* irrFunc .* z_eu;

        nppG1dg(:, :, j) = npp;
    end

    %% Output npp
    filename_ot = ['npp_vgpm_', yyyy];
    save( fullfile( PATH_OT_npp, filename_ot ), 'nppG1dg', 'lat', 'lon', '-v7.3' )

end


%%
function dayL = computeDayLength( lat, r, c, mm, yyyy )
% Computes monthly-mean day length [decimal hours] using the same solar
% declination formulation as computeMonthlyOptics() in
% grid_npp_1M_1dg_toED_ot_modelMLD.m
% lat: Latitude vector [r, 1]
% r, c: Grid dimensions
% mm: Month
% yyyy: Year
    PI = 3.141592653589793;

    num_days = eomday( yyyy, mm );
    nDays = yeardays( yyyy );

    lat_grid = repmat( lat, 1, c );
    dl = nan( r, c, num_days );

    for i = 1 : num_days
        yd = day( datetime( yyyy, mm, i ), 'dayofyear' );

        % Declination (radians)
        decl = 23.5 * cos( 2 * PI * (yd - 172) / nDays ) * PI / 180;

        DL = -1 * tan( lat_grid * PI / 180 ) * tan( decl );
        DL = min( max( DL, -1 ), 1 );
        dl(:, :, i) = acos( DL ) / PI; % Fraction of day
    end

    dayL = nanmean( dl, 3 ) * 24; % Fraction of day -> decimal hours
end
