function script_DEMcharacteristics

dem_path_srtm3 = 'C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_centroids'; %srtm3
dem_path_srtm1 = 'C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_srtm1_centroids'; %srtm1
dem_path_alti3d = 'C:\Users\Simon Rölli\Desktop\data\centroids_large\_LS_Sarnen_interpol_alti3d';%alti3d ca. 7*10m
dem_path_alti3d_2m = 'C:\Users\Simon Rölli\Desktop\data\centroids_large\_LS_Sarnen_alti3d_2m';%alti3d ca. 2x3m
dem_path_alti3d_orig = 'C:\Users\Simon Rölli\Desktop\data\centroids_large\_LS_Sarnen_alti3d_original';%alti3d org 2x2m

% dem_path_srtm3 = 'C:\Users\Simon Rölli\Desktop\data\centroids_large\_LS_Sarnen_srtm3'; %srtm3
% dem_path_srtm1 = 'C:\Users\Simon Rölli\Desktop\data\centroids_large\_LS_Sarnen_srtm1'; %srtm1

%SRTM3
centroids = load(dem_path_srtm3,'centroids');
centroids = centroids.centroids;
n_lon = numel(unique(centroids.lon));
n_lat = numel(unique(centroids.lat));
lon_srtm3 = reshape(centroids.lon,n_lat,n_lon);
lat_srtm3 = reshape(centroids.lat,n_lat,n_lon);
elevation_srtm3 = reshape(centroids.elevation_m,n_lat,n_lon);

%SRTM1
centroids = load(dem_path_srtm1,'centroids');
centroids = centroids.centroids;
n_lon = numel(unique(centroids.lon));
n_lat = numel(unique(centroids.lat));
lon_srtm1 = reshape(centroids.lon,n_lat,n_lon);
lat_srtm1 = reshape(centroids.lat,n_lat,n_lon);
elevation_srtm1 = reshape(centroids.elevation_m,n_lat,n_lon);

% %ALTI3D
% centroids = load(dem_path_alti3d,'centroids');
% centroids = centroids.centroids;
% n_lon = numel(unique(centroids.lon));
% n_lat = numel(unique(centroids.lat));
% lon_alti3d = reshape(centroids.lon,n_lat,n_lon);
% lat_alti3d = reshape(centroids.lat,n_lat,n_lon);
% elevation_alti3d = reshape(centroids.elevation_m,n_lat,n_lon);
% 
% %ALTI3D 2x3m
% centroids = load(dem_path_alti3d_2m,'centroids');
% centroids = centroids.centroids;
% n_lon = numel(unique(centroids.lon));
% n_lat = numel(unique(centroids.lat));
% lon_alti3d_2m = reshape(centroids.lon,n_lat,n_lon);
% lat_alti3d_2m = reshape(centroids.lat,n_lat,n_lon);
% elevation_alti3d_2m = reshape(centroids.elevation_m,n_lat,n_lon);

flowacc = climada_ls_flowacc(lon_srtm1,lat_srtm1,elevation_srtm1,0);
figure 
surface(elevation_srtm1,log(flowacc),'LineStyle','none')

flowacc2 = climada_ls_flowacc(lon_srtm1,lat_srtm1,elevation_srtm1,1);
figure 
surface(elevation_srtm1,log(flowacc2),'LineStyle','none')

%[slope,aspect,area,TWI] = climada_centroids_scores(lon_srtm3,lat_srtm3,dem,TWI,topoTWI)

flowacc = climada_ls_flowacc(mult_flow,0);

DEM = GRIDobj(lon_srtm3,lat_srtm3,elevation_srtm3);
flowacc2 = climada_ls_flowacc(DEM,1);






end