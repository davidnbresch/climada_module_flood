function ls_test2()
%function to test stuff... can be deleted afterwards

%load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_centroids.mat')
load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_srtm1_centroids.mat')

n_lon = numel(unique(centroids.lon));
n_lat = numel(unique(centroids.lat));
lon = reshape(centroids.lon,n_lat,n_lon);
lat = reshape(centroids.lat,n_lat,n_lon);
z = reshape(centroids.elevation_m,n_lat,n_lon);

z = deminpaint(z);


DEM = GRIDobj(lon,lat,flipud(z));
DEM = fillsinks(DEM);
FD = FLOWobj(DEM,'multi');

%DB = drainagebasins(FD,80000);
A = flowacc(FD);
%k = identifyflats(DEM);


figure
surface(z,log(A.Z),'LineStyle','none')


%with old version 1.06
% z = fillsinks(z);
% fl = flowdir(lon,lat,z,'routeflats','route');
% a = flowacc(fl,[n_lat n_lon]);
% figure
% surface(z,log(a),'LineStyle','none')

disp('hier')


end