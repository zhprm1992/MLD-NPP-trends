% Used for: Fig. 3 helper function
% Original: Code\npp_cbpm\plotClimNPP.m (data paths updated to I:\Prj; computation unchanged)
function h = plotClimNPP( ax, MLD, lat, color )
% Plot climatological mld over period of 2003 - 2022
% Haipeng Zhao
% 10/16/2025

MLD_median = nanmedian( MLD, 2 );

% x0 = 1075.4 ; y0 = 50.2; width = 364.8; height = 835.6;
% set( gcf, 'position', [x0, y0, width, height] )
h = plot( ax, MLD_median, lat, 'LineWidth', 1.2, 'Color', color );
% set(gca,'ytick',[], 'FontWeight', 'bold', 'FontSize', 15 )
c = char('S', '', 'N');
yticklabels(arrayfun(@(x) [num2str(abs(x)), char(176), c(sign(x) + 2)], get( ax, 'YTick' ), 'UniformOutput', false));
% xlabel( 'MLD (m)', 'fontweight', 'bold' )

Q1 = prctile(MLD, 10, 2);
Q3 = prctile(MLD, 90, 2);

% ax.FontSize = 15; 
ax.YLim = [-60, 60];
ax.XLim = [0, 1500];

hold on
plot( ax, Q1, lat, 'color', color )
plot( ax, Q3, lat, 'color', color )

% xline( ax, ax.XLim(2) )
% yline( ax, ax.YLim(2) )

% title( 'Spring' )

box off
grid off

end