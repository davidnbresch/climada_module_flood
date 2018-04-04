function write_ascii()
%just to transform lat,lon,elvation matrices in ascii file with meters as
%coordinates (same resolution in x and y direction) --> therefore, interpolation necessary
%can be used in FlowR afterwards --> adjust header lat lon limits to 200000
%600000

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
% cellsizelonlat = lat(1)-lat(2);
% cellsizexyinterp = Y_i(1)-Y_i(2);
% lat2 = lat+abs(cellsizelonlat/2);
% lon2 = lon+abs(cellsizelonlat/2);
% Y_i2 = Y_i+abs(cellsizexyinterp/2);
% X_i2 = X_i+abs(cellsizexyinterp/2);

%%%%%%%%%%%%
%write ascii files

%DEM = GRIDobj(lon,lat,elevation);
%DEM = GRIDobj(X_i,Y_i,elevation_i);
%GRIDobj2ascii(DEM);

%%%%%%%%%%%%%%%%
%write tiff file 

%DEM
%DEM = GRIDobj(lon,lat,elevation)
%GRIDobj2geotiff(DEM);

%source file
source_area = zeros(size(elevation));
source_area(253,51) = 1;
source_area_i = interp2(X,Y,source_area,X_i,Y_i,'nearest');
SCR = GRIDobj(X_i,Y_i,source_area_i);
GRIDobj2ascii(SCR);









%K = GRIDobj(lon,lat,full(k));
%GRIDobj2geotiff(K)

% rasterwrite('source.asc',lat,lon,full(k));
% rasterwrite('elevation.asc',lat,lon,elevation); 

end