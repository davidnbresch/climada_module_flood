function ls_test2()
%function to test stuff... can be deleted afterwards

load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_centroids.mat')
%load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_srtm1_centroids.mat')

n_lon = numel(unique(centroids.lon));
n_lat = numel(unique(centroids.lat));
lon = reshape(centroids.lon,n_lat,n_lon);
lat = reshape(centroids.lat,n_lat,n_lon);
elevation = reshape(centroids.elevation_m,n_lat,n_lon);


DEM = GRIDobj(lon,lat,flipud(elevation));
DEMf = fillsinks(DEM);
FD = FLOWobj(DEMf);
A = flowacc(FD);



disp('hier')


end