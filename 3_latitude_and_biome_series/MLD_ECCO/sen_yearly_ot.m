% Used for: Sen slope of annual MLD-ECCO per grid cell, sens_slope_yearly (Fig. 2c map)
% Original: Code\mld_ECCO\sen_yearly_ot.m (data paths updated to I:\Prj; computation unchanged)
% Calculate sen's slope for annual MLD 
% Haipeng
% 12/07/2023
% 11/13/2024

clear
clc

PATH_IN = 'I:\Prj\MLD_ECCO\biomes_1Y_1dg';
PATH_OT = 'I:\Prj\MLD_ECCO\sens_slope_yearly';

datafile = dir( [PATH_IN, '\', '*.mat'] );
filenum = length( datafile );
for i = 1 : filenum
    %i = 2; % CHECK HERE too
    regfile = datafile(i);
    filename = regfile.name;
    dash_loc = find( filename == '_' );
    biomeName = filename(1 : dash_loc - 1);
    disp( biomeName )
    d = load( [PATH_IN, '\', filename] );
    u = size( d.mld );
    slope = nan( u(1), 1 );
    p_value = nan( u(1), 1 );
    x = 2003 : 2022;
    for j = 1 : u(1)
        y = d.mld(j, :);
        if sum( isnan( y ) ) > 2
            continue
        end
        x = 1 : length(y);
        datain = [x', y'];
        [~, ~, ~, sig, ~, ~, ~, sen, ~, plotofslope, ~] = ktaub(datain, 0.9999, 1);       
%         [~, z_score, sig, ~] = Modified_MannKendall_YW(x', y', 0.05, 0.05);

        if ~isnan(plotofslope)
            % slope(j) = single( (plotofslope(2, 2) - plotofslope(1, 2)) / plotofslope(1, 2) * 100 );
            slope(j) = single( (plotofslope(2, 2) - plotofslope(1, 2)) * 10 );
        end
        p_value(j) = sig;
    end
    idx = d.idx;
    lat = d.lat;
    lon = d.lon;
    save( [PATH_OT, '\', filename], 'slope', 'p_value', 'lon', 'lat', 'idx' )
end