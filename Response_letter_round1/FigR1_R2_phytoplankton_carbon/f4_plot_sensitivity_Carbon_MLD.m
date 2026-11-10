% Used for: round-1 letter Fig. R1 (legend moved by hand)
% Original: Manuscript\MLD\manuscript\nature\NCC\R1\code for revision\sensitivty_carbon_mld\f4_plot_sensitivity_Carbon_MLD.m (data paths updated to I:\Prj; computation unchanged)
% Plot Carbon and MLD sensitivity vs Lat (median version)

clear
clc

load senMLD_carbon_v2.mat
figure('Units','centimeters','Position',[2 2 30 15]);
t = tiledlayout(1, 4);
t.Padding = 'compact';

cmp = slanCM( 'tab20c' );
cmp = unique( cmp, "rows", 'stable' );
c(1, :) = cmp(7, :);
c(2, :) = cmp(5, :);
c(3, :) = cmp(2, :);
c(4, :) = cmp(17, :);

c(1, :) = [120, 120, 120] / 255;
c(2, :) = [230, 159, 0]   / 255;
c(3, :) = [204, 121, 167]   / 255;

senMLD_NEMO(isnan( senMLD_NEMO )) = 0;
senMLD_ECCO(isnan( senMLD_ECCO )) = 0;
senMLD_HYCOM(isnan( senMLD_HYCOM )) = 0;


%% Bar graph (Spring)
nexttile
sen = cat( 2, senMLD_NEMO(:, 1, 1), senMLD_ECCO(:, 1, 1), senMLD_HYCOM(:, 1, 1) );
err = cat( 2, senMLD_NEMO(:, 1, 2), senMLD_ECCO(:, 1, 2), senMLD_HYCOM(:, 1, 2) );
n = cat( 2, senMLD_NEMO(:, 1, 3), senMLD_ECCO(:, 1, 3), senMLD_HYCOM(:, 1, 3) );

y = 1 : 6;
x = flip( sen, 1 );
h = barh( y, x );
drawnow
% ylim( [-3 3] )
for i =1 : 3
    h(i).FaceColor = c(i, :);
    h(i).EdgeColor = c(i, :);
    h(i).LineWidth = 0.8;
    h(i).FaceAlpha = 0.9;
end

hold on;
for k = 1 : size( x, 2 )
    % get x positions per group
    ypos = flip( h(k).XData + h(k).XOffset );
    for t = 1 : 6
        if k == 1
            pt = senPt_nemo{t, 1};
        elseif k == 2
            pt = senPt_ecco{t, 1};
        elseif k == 3
            pt = senPt_hycom{t, 1};
        end
        carbonPt = pt(:, 1) ./ pt(:, 2);
        carbonPt = carbonPt( ~isinf( carbonPt ) );
        if length( carbonPt ) == 1
            yposPt = ypos(t);
        else
            yposPt = linspace( ypos(t) - h(k).BarWidth / 18, ypos(t) + h(k).BarWidth / 18, length( carbonPt ) );
        end
        scatter( carbonPt, yposPt, 20, 'MarkerFaceColor', cmp(k, :), 'MarkerFaceAlpha', 0.1, 'MarkerEdgeColor', c(k, :), 'LineWidth', 1 )
    end
    % draw errorbar
    idx = err(:, k) ~= 0;
    errorbar( sen( idx, k ), ypos(idx)', err(idx, k), 'horizontal', 'LineStyle', 'none', ...
        'Color', 'k', 'LineWidth', 1 );
end

leg = legend( gca, {'NEMO', 'ECCO', 'HYCOM'}, 'Orientation', 'horizontal' );
title( 'Spring' )
set( gca, 'YTickLabel', {'40 - 60S', '20 - 40S', '0 - 20S', '0 - 20N', '20 - 40N', '40 - 60N'} )
% set( gca, 'YTickLabel', {'', '', '', '', '', ''} )
xlim( [-3 3] )


%% Bar graph (SUMMER)
nexttile
sen = cat( 2, senMLD_NEMO(:, 2, 1), senMLD_ECCO(:, 2, 1), senMLD_HYCOM(:, 2, 1) );
err = cat( 2, senMLD_NEMO(:, 2, 2), senMLD_ECCO(:, 2, 2), senMLD_HYCOM(:, 2, 2) );
n = cat( 2, senMLD_NEMO(:, 2, 3), senMLD_ECCO(:, 2, 3), senMLD_HYCOM(:, 2, 3) );

y = 1 : 6;
x = flip( sen, 1 );
h = barh( y, x );
drawnow
% ylim( [-3 3] )
for i =1 : 3
    h(i).FaceColor = c(i, :);
    h(i).EdgeColor = c(i, :);
    h(i).LineWidth = 0.8;
    h(i).FaceAlpha = 0.9;
end

hold on;
for k = 1 : size( x, 2 )
    % get x positions per group
    ypos = flip( h(k).XData + h(k).XOffset );
    for t = 1 : 6
        if k == 1
            pt = senPt_nemo{t, 2};
        elseif k == 2
            pt = senPt_ecco{t, 2};
        elseif k == 3
            pt = senPt_hycom{t, 2};
        end
        carbonPt = pt(:, 1) ./ pt(:, 2);
        carbonPt = carbonPt( ~isinf( carbonPt ) );
        if length( carbonPt ) == 1
            yposPt = ypos(t);
        else
            yposPt = linspace( ypos(t) - h(k).BarWidth / 18, ypos(t) + h(k).BarWidth / 18, length( carbonPt ) );
        end
        scatter( carbonPt, yposPt, 20, 'MarkerFaceColor', cmp(k, :), 'MarkerFaceAlpha', 0.1, 'MarkerEdgeColor', c(k, :), 'LineWidth', 1 )
    end
    % draw errorbar
    idx = err(:, k) ~= 0;
    errorbar( sen( idx, k ), ypos(idx)', err(idx, k), 'horizontal', 'LineStyle', 'none', ...
        'Color', 'k', 'LineWidth', 1 );
end

title( 'Summer' )
% set( gca, 'YTickLabel', {'40 - 60S', '20 - 40S', '0 - 20S', '0 - 20N', '20 - 40N', '40 - 60N'} )
set( gca, 'YTickLabel', '' )

% set( gca, 'YTickLabel', {'', '', '', '', '', ''} )
xlim( [-3 3] )


%% Bar graph (AUTUMN)
nexttile
sen = cat( 2, senMLD_NEMO(:, 3, 1), senMLD_ECCO(:, 3, 1), senMLD_HYCOM(:, 3, 1) );
err = cat( 2, senMLD_NEMO(:, 3, 2), senMLD_ECCO(:, 3, 2), senMLD_HYCOM(:, 3, 2) );
n = cat( 2, senMLD_NEMO(:, 3, 3), senMLD_ECCO(:, 3, 3), senMLD_HYCOM(:, 3, 3) );

y = 1 : 6;
x = flip( sen, 1 );
h = barh( y, x );
drawnow
% ylim( [-3 3] )
for i =1 : 3
    h(i).FaceColor = c(i, :);
    h(i).EdgeColor = c(i, :);
    h(i).LineWidth = 0.8;
    h(i).FaceAlpha = 0.9;
end

hold on;
for k = 1 : size( x, 2 )
    % get x positions per group
    ypos = flip( h(k).XData + h(k).XOffset );
    for t = 1 : 6
        if k == 1
            pt = senPt_nemo{t, 3};
        elseif k == 2
            pt = senPt_ecco{t, 3};
        elseif k == 3
            pt = senPt_hycom{t, 3};
        end
        carbonPt = pt(:, 1) ./ pt(:, 2);
        carbonPt = carbonPt( ~isinf( carbonPt ) );
        if length( carbonPt ) == 1
            yposPt = ypos(t);
        else
            yposPt = linspace( ypos(t) - h(k).BarWidth / 18, ypos(t) + h(k).BarWidth / 18, length( carbonPt ) );
        end
        scatter( carbonPt, yposPt, 20, 'MarkerFaceColor', cmp(k, :), 'MarkerFaceAlpha', 0.1, 'MarkerEdgeColor', c(k, :), 'LineWidth', 1 )
    end
    % draw errorbar
    idx = err(:, k) ~= 0;
    errorbar( sen( idx, k ), ypos(idx)', err(idx, k), 'horizontal', 'LineStyle', 'none', ...
        'Color', 'k', 'LineWidth', 1 );
end

% xlabel( 'Standarized sensitivity', 'FontWeight','bold' )
title( 'Autumn' )
% set( gca, 'YTickLabel', {'40 - 60S', '20 - 40S', '0 - 20S', '0 - 20N', '20 - 40N', '40 - 60N'} )
set( gca, 'YTickLabel', '' )

% set( gca, 'YTickLabel', {'', '', '', '', '', ''} )
xlim( [-3 3] )


%% Bar graph (WINTER)
nexttile
sen = cat( 2, senMLD_NEMO(:, 4, 1), senMLD_ECCO(:, 4, 1), senMLD_HYCOM(:, 4, 1) );
err = cat( 2, senMLD_NEMO(:, 4, 2), senMLD_ECCO(:, 4, 2), senMLD_HYCOM(:, 4, 2) );
n = cat( 2, senMLD_NEMO(:, 4, 3), senMLD_ECCO(:, 4, 3), senMLD_HYCOM(:, 4, 3) );

y = 1 : 6;
x = flip( sen, 1 );
h = barh( y, x );
drawnow
% ylim( [-3 3] )
for i =1 : 3
    h(i).FaceColor = c(i, :);
    h(i).EdgeColor = c(i, :);
    h(i).LineWidth = 0.8;
    h(i).FaceAlpha = 0.9;
end

hold on;
for k = 1 : size( x, 2 )
    % get x positions per group
    ypos = flip( h(k).XData + h(k).XOffset );
    for t = 1 : 6
        if k == 1
            pt = senPt_nemo{t, 4};
        elseif k == 2
            pt = senPt_ecco{t, 4};
        elseif k == 3
            pt = senPt_hycom{t, 4};
        end
        carbonPt = pt(:, 1) ./ pt(:, 2);
        carbonPt = carbonPt( ~isinf( carbonPt ) );
        if length( carbonPt ) == 1
            yposPt = ypos(t);
        else
            yposPt = linspace( ypos(t) - h(k).BarWidth / 18, ypos(t) + h(k).BarWidth / 18, length( carbonPt ) );
        end
        scatter( carbonPt, yposPt, 20, 'MarkerFaceColor', cmp(k, :), 'MarkerFaceAlpha', 0.1, 'MarkerEdgeColor', c(k, :), 'LineWidth', 1 )
    end
    % draw errorbar
    idx = err(:, k) ~= 0;
    errorbar( sen( idx, k ), ypos(idx)', err(idx, k), 'horizontal', 'LineStyle', 'none', ...
        'Color', 'k', 'LineWidth', 1 );
end

title( 'Winter' )
% set( gca, 'YTickLabel', {'40 - 60S', '20 - 40S', '0 - 20S', '0 - 20N', '20 - 40N', '40 - 60N'} )
set( gca, 'YTickLabel', '' )

% set( gca, 'YTickLabel', {'', '', '', '', '', ''} )
xlim( [-3 3] )

