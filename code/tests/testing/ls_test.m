function ls_test()
%function to test stuff... can be deleted afterwards

path = 'C:\Users\Simon Rölli\Desktop\climada\climada_data\alti3d\alti3d.tif';

load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_srtm1_centroids.mat')

[SRTM,srtm_info] = climada_srtm1_get([8.111086 8.341411 46.815369 46.946240],'','','',1);

%get gridded datasets
n_lon = numel(unique(centroids.lon));
n_lat = numel(unique(centroids.lat));
lon = reshape(centroids.lon,n_lat,n_lon);
lat = reshape(centroids.lat,n_lat,n_lon);
elevation = reshape(centroids.elevation_m,n_lat,n_lon);

DEM = GRIDobj(path);

[lon2,lat2] = getcoordinates(DEM);
[lon2,lat2] = meshgrid(lon2,lat2);

deg_km = 111.12;

dy = abs(min(diff(lat(:,1)))*(deg_km * 1000));
dx = abs(min(diff(lon(1,:)))*cosd(mean(lat(:,1)))*(deg_km * 1000)); 
dxdy = sqrt(dx^2+dy^2);

dy2 = abs(min(diff(lat2(:,1)))*(deg_km * 1000));
dx2 = abs(min(diff(lon2(1,:)))*cosd(mean(lat2(:,1)))*(deg_km * 1000)); 
dxdy2 = sqrt(dx2^2+dy2^2);







disp('hier')


end