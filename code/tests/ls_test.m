function ls_test()
%function to test stuff... can be deleted afterwards

% load('C:\Users\Simon Rölli\Desktop\data\z.mat','z');
% load('C:\Users\Simon Rölli\Desktop\data\lon.mat','lon');
% load('C:\Users\Simon Rölli\Desktop\data\lat.mat','lat');

S = shaperead('C:\Users\Simon Rölli\Documents\ArcGIS\Rutschungen_OW\extracted_polygones.shp');

%load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_centroids.mat')
load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_srtm1_centroids.mat')

%load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_hazard.mat')
load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_srtm1_hazard.mat')

n_lon = numel(unique(centroids.lon));
n_lat = numel(unique(centroids.lat));
lon = reshape(centroids.lon,n_lat,n_lon);
lat = reshape(centroids.lat,n_lat,n_lon);
elevation = reshape(centroids.elevation_m,n_lat,n_lon);

latlim = [lat(1,1) lat(n_lat,1)];
lonlim = [lon(1,1) lon(1,n_lon)];

%filled raster from polygones
[Z,R] = vec2mtx([S.Y],[S.X],3600,'filled');
%%%nicht ganz sicher ob eine Zeile verrutscht
[Z,R] = vec2mtx([S.Y],[S.X],3600,latlim,lonlim,'filled');
%[Z,R] = vec2mtx([S.Y],[S.X],zeros(360,360),[3600,latlim(2),lonlim(2)]);

lonZ = [1:numel(Z(1,:))]*(1/R(1))+R(3);
latZ = R(2)-[1:numel(Z(:,1))]*(1/R(1));

[lonZ,latZ] = meshgrid(lonZ,latZ);
latZ = flipud(latZ);

DEM = GRIDobj(lonZ,latZ,Z);
GRIDobj2geotiff(DEM);



for i = 1:numel(S)
    %[Z, R] = vec2mtx(S(i).Y,S(i).X,3600,, lonlim)
end

geotif

disp('hier')


end