function dem4FlowR()
%not a climada function
%can be used to transform dem with lat lon in x and y [m]
%coordinates...because flowR ask for the same resolution in x and y
%direction, after the lon/lat is transformed into x/y an interpolation is
%conducted to get the same resolution.
%the new dem should then be used to compare individual slides with FlowR
%(climada needs lat/lon --> calculate back such that it uses the same DEM
%for the comparison)

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

%derive raster grid in meters, x coordinates need to be adjsted to have
%same resolution as y coordinates --> therefore DEM interpolation from
%original x to adjusted x need to be conducted
deg_km = 111.12;
dy = abs(min(diff(lat(:,1)))*(deg_km * 1000));
dx = abs(min(diff(lon(1,:)))*cosd(mean(lat(:,1)))*(deg_km * 1000));
x = [1:360]*dx;
nX = floor(x(numel(x))/dy);
x_i = [1:nX]*dy;
y = ([1:360]*dy)';
[X,Y] = meshgrid(x,y);
[X_i,Y_i] = meshgrid(x_i,y);
elevation_i = interp2(X,Y,elevation,X_i,Y_i,'linear');

%I think coordinates are already pointing to middle of gridcell, therefore
%add a half length of cellsize in lat/lon direction. Because GRIDobj2ascii
%substract same amount afterwards
cellsizelonlat = lat(1)-lat(2);
cellsizexyinterp = Y_i(1)-Y_i(2);
lat2 = lat+abs(cellsizelonlat/2);
lon2 = lon+abs(cellsizelonlat/2);
Y_i2 = Y_i+abs(cellsizexyinterp/2);
X_i2 = X_i+abs(cellsizexyinterp/2);

%%%%%%%%%%%%
%write ascii files

DEM = GRIDobj(lon2,lat2,elevation);
%DEM = GRIDobj(X,Y,elevation); %not working because not same lat/lon resolution
%DEM = GRIDobj(X_i2,Y_i2,elevation_i);
%GRIDobj2ascii(DEM);

%%%%%%%%%%%%%%%%
%write tiff file 

%DEM
DEM = GRIDobj(lon,lat,elevation)
GRIDobj2geotiff(DEM);

%source file
source_area = reshape(hazard.intensity(1,:),n_lat,n_lon);








%K = GRIDobj(lon,lat,full(k));
%GRIDobj2geotiff(K)

% rasterwrite('source.asc',lat,lon,full(k));
% rasterwrite('elevation.asc',lat,lon,elevation); 

end