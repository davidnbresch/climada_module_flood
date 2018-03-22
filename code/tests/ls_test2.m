function ls_test2()
%function to test stuff... can be deleted afterwards

%load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_centroids.mat')
load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_srtm1_centroids.mat')

%load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_hazard.mat')
load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_srtm1_hazard.mat')

n_lon = numel(unique(centroids.lon));
n_lat = numel(unique(centroids.lat));
lon = reshape(centroids.lon,n_lat,n_lon);
lat = reshape(centroids.lat,n_lat,n_lon);
elevation = reshape(centroids.elevation_m,n_lat,n_lon);

elevation = deminpaint(elevation);
elevation = fillsinks(elevation);

%test rasterwrite
deg_km = 111.12;
dy = abs(min(diff(lat(:,1)))*(deg_km * 1000));
dx = abs(min(diff(lon(1,:)))*cosd(mean(lat(:,1)))*(deg_km * 1000));
x = [1:360]*dx;
nX = floor(x(numel(x))/dy);
x_i = [1:nX]*dy;
y = ([1:360]*dy)';
[X,Y] = meshgrid(x,y);
[X_i,Y_i] = meshgrid(x_i,y);
% rasterwrite('asciidemlatlon.asc',lat,lon,elevation);
% rasterwrite('asciidemYX.asc',Y,X,elevation);

elevation_i = interp2(X,Y,elevation,X_i,Y_i,'linear');

source_area = reshape(hazard.intensity(1,:),n_lat,n_lon);

DEM = GRIDobj(lon,lat,elevation);
%GRIDobj2geotiff(DEM)
GRIDobj2ascii(DEM);

K = GRIDobj(lon,lat,full(k));
%GRIDobj2geotiff(K)

% rasterwrite('source.asc',lat,lon,full(k));
% rasterwrite('elevation.asc',lat,lon,elevation); 

end