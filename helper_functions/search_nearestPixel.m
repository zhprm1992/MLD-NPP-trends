% Used for: nearest 1/12 or LLC270 cell to each 1 deg grid point (regridding index scripts in 1_MLD_computation)
% Original: Code\aerosol\search_nearestPixel.m (data paths updated to I:\Prj; computation unchanged)
function idx = search_nearestPixel( lon_station, lat_station, lon_pixels, lat_pixels, boundary )
% Search for the nearest pixel to fixed station
%
%  search_nearestPixel( lon_station, lat_station, lon_pixels, lat_pixels, 0.06 )
%
% 9/9/2022
% 9/6/2023
% Haipeng Zhao

idx = [1 : length( lon_pixels )]';
idx_lat = lat_pixels > lat_station - boundary & lat_pixels < lat_station +boundary;
idx_lon = lon_pixels > lon_station - boundary & lon_pixels < lon_station + boundary;
idx1 = idx_lat & idx_lon;

lat = lat_pixels(idx1);
lon = lon_pixels(idx1);
idx = idx(idx1);

pt1 = [lat_station, lon_station];
dist = zeros( length( lat ), 1 );
for k = 1 : length( dist )
    pt2 = [lat(k), lon(k)];
    dist(k) = deg2km( distance( pt1, pt2 ) );
end
[d, idx2] = min( dist );
idx = idx(idx2);

 if isempty( lat | lon )
    idx = nan;
 end       

end