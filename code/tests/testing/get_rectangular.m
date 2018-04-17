function [elevation2,x2,y2,lon2,lat2] = get_rectangular(lon,lat,elevation)
%just to transform lat,lon,elvation matrices in ascii file with meters as
%coordinates (same resolution in x and y direction) --> therefore, interpolation necessary
%can be used in FlowR afterwards --> adjust header lat lon limits to 200000
%600000

% %load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_centroids.mat')
% load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_srtm1_centroids.mat')
% 
% %load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_hazard.mat')
% load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_srtm1_hazard.mat')
% 
% n_lon = numel(unique(centroids.lon));
% n_lat = numel(unique(centroids.lat));
% lon = reshape(centroids.lon,n_lat,n_lon);
% lat = reshape(centroids.lat,n_lat,n_lon);
% elevation = reshape(centroids.elevation_m,n_lat,n_lon);

%elevation = deminpaint(elevation);
elevation = fillsinks(elevation);

%derive raster grid in meters, x coordinates need to be adjsted to have
%same resolution as y coordinates --> therefore DEM interpolation from
%original x to adjusted x need to be conducted
deg_km = 111.12;
dlat = abs(min(diff(lat(:,1))));
dlon = abs(min(diff(lon(1,:))));
dy = dlat*(deg_km * 1000);
dx = dlon*cosd(mean(lat(:,1)))*(deg_km * 1000);

%x y vector of original file
x = [1:numel(lon(1,:))]*dx;
y = ([1:numel(lat(:,1))]*dy)';
[x,y] = meshgrid(x,y);

%resolution in degree to have same metric resolution dy also in x direction
dlat2 = dlat;
dlon2 = dy/(cosd(mean(lat(:,1)))*(deg_km * 1000));
dy2 = dy;
dx2 = dlon2*cosd(mean(lat(:,1)))*(deg_km * 1000);

%calculate how many centroids we need in x direction
nx2 = floor(x(1,numel(x(1,:)))/dx2);

%get new coordinates
%x and y coordinates
y2 = y(:,1);
x2 = [1:nx2]*dx2;
[x2,y2] = meshgrid(x2,y2);
%lon lat coordinates
lat2 = lat(:,1);
lon2 = lon(1,1) + [1:nx2]*dlon2;
[lon2,lat2] = meshgrid(lon2,lat2);

%interpolate elevation data
elevation2 = interp2(lon,lat,elevation,lon2,lat2,'linear');


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
% source_area = zeros(size(elevation));
% source_area(253,51) = 1;
% source_area_i = interp2(X,Y,source_area,X_i,Y_i,'nearest');
% SCR = GRIDobj(X_i,Y_i,source_area_i);
% GRIDobj2ascii(SCR);









%K = GRIDobj(lon,lat,full(k));
%GRIDobj2geotiff(K)

% rasterwrite('source.asc',lat,lon,full(k));
% rasterwrite('elevation.asc',lat,lon,elevation); 

end