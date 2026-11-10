% Used for: trends per latitude row and month, NPPtrends.mat (Fig. 4) and Fig. S6; the save line at the end is switched on
% Original: Code\npp_cbpm\fs6_heatmap_monthlyNPPtrends_ot.m (data paths updated to I:\Prj; computation unchanged)
% lot heat-map latitudinal trends of chl and four drivers (60s-60n) vs Month, 
% Haipeng Zhao
% 2/23/2024
% 1/27/2026

clear
clc

path = 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\nppLat_monthly_nemo.mat';
d = load( path );
idx = d.lat < 60 & d.lat > -60;
lat = d.lat(idx);
nppLat_nemo = d.nppLat(idx, :);

path = 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\nppLat_monthly_ecco.mat';
d = load( path );
idx = d.lat < 60 & d.lat > -60;
nppLat_ecco = d.nppLat(idx, :);

path = 'I:\Prj\npp_cbpm_1M_1dg_zno3mMU\nppLat_monthly_hycom.mat';
d = load( path );
idx = d.lat < 60 & d.lat > -60;
nppLat_hycom = d.nppLat(idx, :);

path = 'I:\Prj\MLD_1M_NEMO\mldLat_monthly_nemo.mat';
d = load( path );
idx = d.lat < 60 & d.lat > -60;
mldLat_nemo = d.mldLat(idx, :);

path = 'I:\Prj\MLD_ECCO\mldLat_monthly_ecco.mat';
d = load( path );
idx = d.lat < 60 & d.lat > -60;
mldLat_ecco = d.mldLat(idx, :);

path = 'I:\Prj\MLD_HYCOM\mldLat_monthly_hycom.mat';
d = load( path );
idx = d.lat < 60 & d.lat > -60;
mldLat_hycom = d.mldLat(idx, :);

%%
npp_trends_pct_nemo = nan( 12, length( lat ) );
npp_trends_nemo = nan( 12, length( lat ) );
npp_plotofslope_nemo = nan( 12, length( lat ) );
npp_trends_sig_nemo = nan( 12, length( lat ) );

npp_trends_pct_ecco = nan( 12, length( lat ) );
npp_trends_ecco = nan( 12, length( lat ) );
npp_plotofslope_ecco = nan( 12, length( lat ) );
npp_trends_sig_ecco = nan( 12, length( lat ) );

npp_trends_pct_hycom = nan( 12, length( lat ) );
npp_trends_hycom = nan( 12, length( lat ) );
npp_plotofslope_hycom = nan( 12, length( lat ) );
npp_trends_sig_hycom = nan( 12, length( lat ) );

mld_trends_pct_nemo = nan( 12, length( lat ) );
mld_trends_nemo = nan( 12, length( lat ) );
mld_plotofslope_nemo = nan( 12, length( lat ) );
mld_trends_sig_nemo = nan( 12, length( lat ) );

mld_trends_pct_ecco = nan( 12, length( lat ) );
mld_trends_ecco = nan( 12, length( lat ) );
mld_plotofslope_ecco = nan( 12, length( lat ) );
mld_trends_sig_ecco = nan( 12, length( lat ) );

mld_trends_pct_hycom = nan( 12, length( lat ) );
mld_trends_hycom = nan( 12, length( lat ) );
mld_plotofslope_hycom = nan( 12, length( lat ) );
mld_trends_sig_hycom = nan( 12, length( lat ) );

fig = figure('Units', 'inches', 'Position', [1, 1, 7, 7.2]);  % [x, y, width, height]

t = tiledlayout(2, 3);
t.TileSpacing = 'loose';
t.Padding = 'compact';

%% npp_nemo
nexttile
for i = 1 : length( lat )
    for j = 1 : 12
        npp_lat = nppLat_nemo(i, d.mm == j);
        if sum( isnan( npp_lat ) ) > 0
            continue
        end
        period_length = size( npp_lat, 2 );
        datain = [(1 : period_length)', npp_lat'];
        [~, ~, r] = regress( datain(:, 2), [(ones( length( datain(:, 1) ), 1 )), datain(:, 1)]  );
        [p_ac_npp(i), ~] = dwtest( r, datain(:, 1) );
        [~, ~, ~, sig, ~, ~, ~, sen, ~, plotofslope, ~] = ktaub(datain, 0.9999, 1);       
        p_value = sig;
        if ~isnan(plotofslope)
            if plotofslope(1, 2) == 0
                npp_trends_pct_nemo(j, i) = 0;
                npp_trends_nemo(j, i) = 0;
                npp_plotofslope_nemo(j, i) = 0;
                npp_trends_sig_nemo(j, i) = p_value;
            else
                npp_trends_pct_nemo(j, i) = single( (plotofslope(2, 2) - plotofslope(1, 2)) / plotofslope(1, 2) * 100 );
                npp_trends_nemo(j, i) = sen;
                npp_plotofslope_nemo(j, i) = plotofslope(1, 2);
                npp_trends_sig_nemo(j, i) = p_value;
            end
        end
    end
end

%
x_label = {'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'};
y = (59 : -1 : -60);
c = char('S', '', 'N');
y_label = arrayfun(@(x) [num2str(abs(x)), char(176), c(sign(x) + 2)], y, 'UniformOutput', false);

cmp = crameri('cork') ;

npp_trends_nemo = rot90( npp_trends_nemo );
npp_trends_pct_nemo = rot90( npp_trends_pct_nemo );
npp_plotofslope_nemo = rot90( npp_plotofslope_nemo );
heatmap( npp_trends_pct_nemo, x_label, y_label, [], [],...
    'MinColorValue', -2, 'MaxColorValue', 2,...
    'ShowAllTicks', 0, 'NaNColor', [0.1 0.1 0.1],...
    'Colormap', cmp, ...
    'TickFontSize', 10, 'UseFigureColormap', 0)
cb = colorbar( 'northoutside' );
cb.Label.String = 'NPP-NEMO trends (% yr^{-1})';
cb.Label.FontWeight = 'bold';

% colorbar
% xlabel( 'Month', 'FontWeight', 'bold' )
% set( gca, 'fontweight', 'bold' )

npp_trends_sig_nemo = rot90( npp_trends_sig_nemo );
x = 1 : 12; y = (1 : 120);
a = 0.5;
for j = 1 : size( npp_trends_sig_nemo, 1 )
    for i = 1 : size( npp_trends_sig_nemo, 2 )
        if npp_trends_sig_nemo(j, i) < 0.001
            text( x(i), y(j) + a, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        elseif npp_trends_sig_nemo(j, i) < 0.01
            text( x(i), y(j) + a, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        elseif npp_trends_sig_nemo(j, i) < 0.05
            text( x(i), y(j) + a, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        end
    end
end

%% npp-ecco
nexttile
for i = 1 : length( lat )
    for j = 1 : 12
        npp_lat = nppLat_ecco(i, d.mm == j);
        if sum( isnan( npp_lat ) ) > 0
            continue
        end
        period_length = size( npp_lat, 2 );
        datain = [(1 : period_length)', npp_lat'];
        [~, ~, r] = regress( datain(:, 2), [(ones( length( datain(:, 1) ), 1 )), datain(:, 1)]  );
        [p_ac_npp(i), ~] = dwtest( r, datain(:, 1) );
        [~, ~, ~, sig, ~, ~, ~, sen, ~, plotofslope, ~] = ktaub(datain, 0.9999, 1);       
        p_value = sig;
        if ~isnan(plotofslope)
            if plotofslope(1, 2) == 0
                npp_trends_pct_ecco(j, i) = 0;
                npp_trends_ecco(j, i) = 0;
                npp_plotofslope_ecco(j, i) = 0;
                npp_trends_sig_ecco(j, i) = p_value;
            else
                npp_trends_pct_ecco(j, i) = single( (plotofslope(2, 2) - plotofslope(1, 2)) / plotofslope(1, 2) * 100 );
                npp_trends_ecco(j, i) = sen;
                npp_plotofslope_ecco(j, i) = plotofslope(1, 2);
                npp_trends_sig_ecco(j, i) = p_value;
            end
        end
    end
end

%
x_label = {'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'};
y = (59 : -1 : -60);
c = char('S', '', 'N');
y_label = arrayfun(@(x) [num2str(abs(x)), char(176), c(sign(x) + 2)], y, 'UniformOutput', false);

cmp = crameri('cork') ;

npp_trends_ecco = rot90( npp_trends_ecco );
npp_trends_pct_ecco = rot90( npp_trends_pct_ecco );
npp_plotofslope_ecco = rot90( npp_plotofslope_ecco );
heatmap( npp_trends_pct_ecco, x_label, [], [], [],...
    'MinColorValue', -2, 'MaxColorValue', 2,...
    'ShowAllTicks', 0, 'NaNColor', [0.1 0.1 0.1],...
    'Colormap', cmp, ...
    'TickFontSize', 10, 'UseFigureColormap', 0)
cb = colorbar( 'northoutside' );
cb.Label.String = 'NPP-ECCO trends (% yr^{-1})';
cb.Label.FontWeight = 'bold';

% colorbar
% xlabel( 'Month', 'FontWeight', 'bold' )
% set( gca, 'fontweight', 'bold' )

npp_trends_sig_ecco = rot90( npp_trends_sig_ecco );
x = 1 : 12; y = (1 : 120);
a = 0.5;
for j = 1 : size( npp_trends_sig_ecco, 1 )
    for i = 1 : size( npp_trends_sig_ecco, 2 )
        if npp_trends_sig_ecco(j, i) < 0.001
            text( x(i), y(j) + a, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        elseif npp_trends_sig_ecco(j, i) < 0.01
            text( x(i), y(j) + a, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        elseif npp_trends_sig_ecco(j, i) < 0.05
            text( x(i), y(j) + a, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        end
    end
end


%% npp-hycom
nexttile
for i = 1 : length( lat )
    for j = 1 : 12
        npp_lat = nppLat_hycom(i, d.mm == j);
        if sum( isnan( npp_lat ) ) > 0
            continue
        end
        period_length = size( npp_lat, 2 );
        datain = [(1 : period_length)', npp_lat'];
        [~, ~, r] = regress( datain(:, 2), [(ones( length( datain(:, 1) ), 1 )), datain(:, 1)]  );
        [p_ac_npp(i), ~] = dwtest( r, datain(:, 1) );
        [~, ~, ~, sig, ~, ~, ~, sen, ~, plotofslope, ~] = ktaub(datain, 0.9999, 1);       
        p_value = sig;
        if ~isnan(plotofslope)
            if plotofslope(1, 2) == 0
                npp_trends_pct_hycom(j, i) = 0;
                npp_trends_hycom(j, i) = 0;
                npp_plotofslope_hycom(j, i) = 0;
                npp_trends_sig_hycom(j, i) = p_value;
            else
                npp_trends_pct_hycom(j, i) = single( (plotofslope(2, 2) - plotofslope(1, 2)) / plotofslope(1, 2) * 100 );
                npp_trends_hycom(j, i) = sen;
                npp_plotofslope_hycom(j, i) = plotofslope(1, 2);
                npp_trends_sig_hycom(j, i) = p_value;
            end
        end
    end
end

%
x_label = {'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'};
y = (59 : -1 : -60);
c = char('S', '', 'N');
y_label = arrayfun(@(x) [num2str(abs(x)), char(176), c(sign(x) + 2)], y, 'UniformOutput', false);

cmp = crameri('cork') ;

npp_trends_hycom = rot90( npp_trends_hycom );
npp_trends_pct_hycom = rot90( npp_trends_pct_hycom );
npp_plotofslope_hycom = rot90( npp_plotofslope_hycom );
heatmap( npp_trends_pct_hycom, x_label, [], [], [],...
    'MinColorValue', -2, 'MaxColorValue', 2,...
    'ShowAllTicks', 0, 'NaNColor', [0.1 0.1 0.1],...
    'Colormap', cmp, ...
    'TickFontSize', 10, 'UseFigureColormap', 0)
set(gca, 'CLim', [-2 2]);
colormap( gca, cmp )
cb = colorbar( 'northoutside' );
cb.Label.String = 'NPP-HYCOM trends (% yr^{-1})';
cb.Label.FontWeight = 'bold';

% colorbar
% xlabel( 'Month', 'FontWeight', 'bold' )
% set( gca, 'fontweight', 'bold' )

npp_trends_sig_hycom = rot90( npp_trends_sig_hycom );
x = 1 : 12; y = (1 : 120);
a = 0.5;
for j = 1 : size( npp_trends_sig_hycom, 1 )
    for i = 1 : size( npp_trends_sig_hycom, 2 )
        if npp_trends_sig_hycom(j, i) < 0.001
            text( x(i), y(j) + a, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        elseif npp_trends_sig_hycom(j, i) < 0.01
            text( x(i), y(j) + a, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        elseif npp_trends_sig_hycom(j, i) < 0.05
            text( x(i), y(j) + a, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        end
    end
end

%% mld-nemo
nexttile
for i = 1 : length( lat )
    for j = 1 : 12
        mld_lat = mldLat_nemo(i, d.mm == j);
        if sum( isnan( mld_lat ) ) > 0
            continue
        end
        period_length = size( mld_lat, 2 );
        datain = [(1 : period_length)', mld_lat'];
        [~, ~, r] = regress( datain(:, 2), [(ones( length( datain(:, 1) ), 1 )), datain(:, 1)]  );
        [p_ac_mld(i), ~] = dwtest( r, datain(:, 1) );
        [~, ~, ~, sig, ~, ~, ~, sen, ~, plotofslope, ~] = ktaub(datain, 0.9999, 1);       
        p_value = sig;
        if ~isnan(plotofslope)
            if plotofslope(1, 2) == 0
                mld_trends_pct_nemo(j, i) = 0;
                mld_trends_nemo(j, i) = 0;
                mld_plotofslope_nemo(j, i) = 0;
                mld_trends_sig_nemo(j, i) = p_value;
            else
                mld_trends_pct_nemo(j, i) = single( (plotofslope(2, 2) - plotofslope(1, 2)) / plotofslope(1, 2) * 100 );
                mld_trends_nemo(j, i) = sen;
                mld_plotofslope_nemo(j, i) = plotofslope(1, 2);
                mld_trends_sig_nemo(j, i) = p_value;
            end
        end
    end
end

%
x_label = {'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'};
y = (59 : -1 : -60);
c = char('S', '', 'N');
y_label = arrayfun(@(x) [num2str(abs(x)), char(176), c(sign(x) + 2)], y, 'UniformOutput', false);

cmp = crameri('roma') ;

mld_trends_nemo = rot90( mld_trends_nemo );
mld_trends_pct_nemo = rot90( mld_trends_pct_nemo );
mld_plotofslope_nemo = rot90( mld_plotofslope_nemo );
heatmap( mld_trends_pct_nemo, x_label, y_label, [], [],...
    'MinColorValue', -2, 'MaxColorValue', 2,...
    'ShowAllTicks', 0, 'NaNColor', [0.1 0.1 0.1],...
    'Colormap', cmp, ...
    'TickFontSize', 10, 'UseFigureColormap', 0)
cb = colorbar( 'northoutside' );
cb.Label.String = 'MLD-NEMO trends (% yr^{-1})';
cb.Label.FontWeight = 'bold';

% colorbar
% xlabel( 'Month', 'FontWeight', 'bold' )
% set( gca, 'fontweight', 'bold' )

mld_trends_sig_nemo = rot90( mld_trends_sig_nemo );
x = 1 : 12; y = (1 : 120);
a = 0.5;
for j = 1 : size( mld_trends_sig_nemo, 1 )
    for i = 1 : size( mld_trends_sig_nemo, 2 )
        if mld_trends_sig_nemo(j, i) < 0.001
            text( x(i), y(j) + a, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        elseif mld_trends_sig_nemo(j, i) < 0.01
            text( x(i), y(j) + a, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        elseif mld_trends_sig_nemo(j, i) < 0.05
            text( x(i), y(j) + a, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        end
    end
end

%% mld-ecco
nexttile
for i = 1 : length( lat )
    for j = 1 : 12
        mld_lat = mldLat_ecco(i, d.mm == j);
        if sum( isnan( mld_lat ) ) > 0
            continue
        end
        period_length = size( mld_lat, 2 );
        datain = [(1 : period_length)', mld_lat'];
        [~, ~, r] = regress( datain(:, 2), [(ones( length( datain(:, 1) ), 1 )), datain(:, 1)]  );
        [p_ac_mld(i), ~] = dwtest( r, datain(:, 1) );
        [~, ~, ~, sig, ~, ~, ~, sen, ~, plotofslope, ~] = ktaub(datain, 0.9999, 1);       
        p_value = sig;
        if ~isnan(plotofslope)
            if plotofslope(1, 2) == 0
                mld_trends_pct_ecco(j, i) = 0;
                mld_trends_ecco(j, i) = 0;
                mld_plotofslope_ecco(j, i) = 0;
                mld_trends_sig_ecco(j, i) = p_value;
            else
                mld_trends_pct_ecco(j, i) = single( (plotofslope(2, 2) - plotofslope(1, 2)) / plotofslope(1, 2) * 100 );
                mld_trends_ecco(j, i) = sen;
                mld_plotofslope_ecco(j, i) = plotofslope(1, 2);
                mld_trends_sig_ecco(j, i) = p_value;
            end
        end
    end
end

%
x_label = {'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'};
y = (59 : -1 : -60);
c = char('S', '', 'N');
y_label = arrayfun(@(x) [num2str(abs(x)), char(176), c(sign(x) + 2)], y, 'UniformOutput', false);

cmp = crameri('roma') ;

mld_trends_ecco = rot90( mld_trends_ecco );
mld_trends_pct_ecco = rot90( mld_trends_pct_ecco );
mld_plotofslope_ecco = rot90( mld_plotofslope_ecco );
heatmap( mld_trends_pct_ecco, x_label, y_label, [], [],...
    'MinColorValue', -2, 'MaxColorValue', 2,...
    'ShowAllTicks', 0, 'NaNColor', [0.1 0.1 0.1],...
    'Colormap', cmp, ...
    'TickFontSize', 10, 'UseFigureColormap', 0)
cb = colorbar( 'northoutside' );
cb.Label.String = 'MLD-ECCO trends (% yr^{-1})';
cb.Label.FontWeight = 'bold';

% colorbar
% xlabel( 'Month', 'FontWeight', 'bold' )
% set( gca, 'fontweight', 'bold' )

mld_trends_sig_ecco = rot90( mld_trends_sig_ecco );
x = 1 : 12; y = (1 : 120);
a = 0.5;
for j = 1 : size( mld_trends_sig_ecco, 1 )
    for i = 1 : size( mld_trends_sig_ecco, 2 )
        if mld_trends_sig_ecco(j, i) < 0.001
            text( x(i), y(j) + a, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        elseif mld_trends_sig_ecco(j, i) < 0.01
            text( x(i), y(j) + a, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        elseif mld_trends_sig_ecco(j, i) < 0.05
            text( x(i), y(j) + a, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        end
    end
end

%% mld-hycom
nexttile
for i = 1 : length( lat )
    for j = 1 : 12
        mld_lat = mldLat_hycom(i, d.mm == j);
        if sum( isnan( mld_lat ) ) > 0
            continue
        end
        period_length = size( mld_lat, 2 );
        datain = [(1 : period_length)', mld_lat'];
        [~, ~, r] = regress( datain(:, 2), [(ones( length( datain(:, 1) ), 1 )), datain(:, 1)]  );
        [p_ac_mld(i), ~] = dwtest( r, datain(:, 1) );
        [~, ~, ~, sig, ~, ~, ~, sen, ~, plotofslope, ~] = ktaub(datain, 0.9999, 1);       
        p_value = sig;
        if ~isnan(plotofslope)
            if plotofslope(1, 2) == 0
                mld_trends_pct_hycom(j, i) = 0;
                mld_trends_hycom(j, i) = 0;
                mld_plotofslope_hycom(j, i) = 0;
                mld_trends_sig_hycom(j, i) = p_value;
            else
                mld_trends_pct_hycom(j, i) = single( (plotofslope(2, 2) - plotofslope(1, 2)) / plotofslope(1, 2) * 100 );
                mld_trends_hycom(j, i) = sen;
                mld_plotofslope_hycom(j, i) = plotofslope(1, 2);
                mld_trends_sig_hycom(j, i) = p_value;
            end
        end
    end
end

%
x_label = {'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'};
y = (59 : -1 : -60);
c = char('S', '', 'N');
y_label = arrayfun(@(x) [num2str(abs(x)), char(176), c(sign(x) + 2)], y, 'UniformOutput', false);

cmp = crameri('roma') ;

mld_trends_hycom = rot90( mld_trends_hycom );
mld_trends_pct_hycom = rot90( mld_trends_pct_hycom );
mld_plotofslope_hycom = rot90( mld_plotofslope_hycom );
heatmap( mld_trends_pct_hycom, x_label, y_label, [], [],...
    'MinColorValue', -2, 'MaxColorValue', 2,...
    'ShowAllTicks', 0, 'NaNColor', [0.1 0.1 0.1],...
    'Colormap', cmp, ...
    'TickFontSize', 10, 'UseFigureColormap', 0)
cb = colorbar( 'northoutside' );
cb.Label.String = 'MLD-HYCOM trends (% yr^{-1})';
cb.Label.FontWeight = 'bold';

% colorbar
% xlabel( 'Month', 'FontWeight', 'bold' )
% set( gca, 'fontweight', 'bold' )

mld_trends_sig_hycom = rot90( mld_trends_sig_hycom );
x = 1 : 12; y = (1 : 120);
a = 0.5;
for j = 1 : size( mld_trends_sig_hycom, 1 )
    for i = 1 : size( mld_trends_sig_hycom, 2 )
        if mld_trends_sig_hycom(j, i) < 0.001
            text( x(i), y(j) + a, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        elseif mld_trends_sig_hycom(j, i) < 0.01
            text( x(i), y(j) + a, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        elseif mld_trends_sig_hycom(j, i) < 0.05
            text( x(i), y(j) + a, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        end
    end
end

save( 'NPPtrends.mat', 'npp_trends_pct_nemo', 'npp_trends_pct_ecco', 'npp_trends_pct_hycom',...
     'mld_trends_pct_nemo', 'mld_trends_pct_ecco', 'mld_trends_pct_hycom',...
      'npp_trends_sig_nemo', 'npp_trends_sig_ecco', 'npp_trends_sig_hycom',...
      'mld_trends_sig_nemo', 'mld_trends_sig_ecco', 'mld_trends_sig_hycom', 'lat' )
