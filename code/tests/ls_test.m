function ls_test()
%function to test stuff... can be deleted afterwards

% load('C:\Users\Simon Rölli\Desktop\data\z.mat','z');
% load('C:\Users\Simon Rölli\Desktop\data\lon.mat','lon');
% load('C:\Users\Simon Rölli\Desktop\data\lat.mat','lat');

S = shaperead('C:\Users\Simon Rölli\Documents\ArcGIS\Rutschungen_OW\selected_forMatlab.shp');

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

%resolution
resol = 3600;
buffer = 1/7200;

lon_72 = [lonlim(1):1/resol:lonlim(2)];
lat_72 = [latlim(1):1/resol:latlim(2)];

[lon_72,lat_72] = meshgrid(lon_72,lat_72);

DEM = GRIDobj(lon_72,flipud(lat_72),zeros(size(lon_72)));
P = polygon2GRIDobj(DEM,S,'TARGET_FID');
% mapshow(lon_72,lat_72,zeros(size(P.Z)),'CData',flipud(P.Z),'DisplayType','surface');
% mapshow(S);



%widening of polygons with 
Sb = S;
for i = 1:numel(S)
    [latb,lonb] = bufferm(S(i).Y,S(i).X,buffer);
    Sb(i).X = lonb;
    Sb(i).Y = latb;
end

Pb = polygon2GRIDobj(DEM,Sb,'TARGET_FID');

merge = P.Z+Pb.Z;

figure
mapshow(lon_72,lat_72,zeros(size(merge)),'CData',flipud(merge),'DisplayType','surface');
mapshow(S);

%use of topotoolbox --> result not good enough --> when interpolating back

[Z,R] = vec2mtx([S.Y],[S.X],resol,latlim,lonlim,'filled');

lonZ = [1:numel(Z(1,:))]*(1/R(1))+R(3)-1/(resol*2);
latZ = R(2)-[1:numel(Z(:,1))]*(1/R(1))+1/(resol*2);

[lonZ,latZ] = meshgrid(lonZ,latZ);
latZ = flipud(latZ);

figure 
mapshow(lonZ,latZ,zeros(size(Z)),'CData',Z,'DisplayType','surface')
mapshow(S)
figure
mapshow(zeros(size(Z)),R,'CData',Z,'DisplayType','surface')
mapshow(S)

interpolat_grid = interp2(lonZ,latZ,Z,lon,lat,'nearest');

figure
mapshow(lon,lat,zeros(size(interpolat_grid)),'CData',interpolat_grid,'DisplayType','surface')
mapshow(S)

%DEM = GRIDobj(lonZ,latZ,Z);
%GRIDobj2geotiff(DEM);

% surface(lon,lat,elevation)
% figure

% for i = 1:numel(S)
%     %[Z, R] = vec2mtx(S(i).Y,S(i).X,3600,, lonlim)
% end
% 
% geotif

disp('hier')


end