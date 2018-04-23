function climada_ls_FPcalibration(lon,lat,S,field,buffer)

% Script to calibrate flow path parameters
% MODULE:
%   flood
% NAME:
%   climada_ls_FPcalibration
% PURPOSE:
%   
% CALLING SEQUENCE:
%   
% EXAMPLE:
%   
% INPUTS: 
%     
% OPTIONAL INPUT PARAMETERS:
%    
% OUTPUTS:
%    
% MODIFICATION HISTORY:
% Thomas Rölli, thomasroelli@gmail.com, 20180417, init

global climada_global
if ~climada_init_vars, return; end

%read in shape file of past shallow landslides
S = shaperead('C:\Users\Simon Rölli\Desktop\data\inventory\slides_forMatlab.shp');

%read in DEM (original alti3d with 2x2m resolution is best choice --> best results when assessing source
%area lenght and width of slides
load('C:\Users\Simon Rölli\Desktop\data\centroids_large\_LS_Sarnen_alti3d_2m.mat','centroids');

n_lon = numel(unique(centroids.lon));
n_lat = numel(unique(centroids.lat));
lon = reshape(centroids.lon,n_lat,n_lon);
lat = reshape(centroids.lat,n_lat,n_lon);
elevation = reshape(centroids.elevation_m,n_lat,n_lon);

%[slope,~,area,~] = climada_centroids_scores(lon,lat,elevation,0);

DEM = GRIDobj(lon,lat,elevation);

%transform polygon in raster
polyraster = climada_ls_poly2raster(lon,lat,S,'OBJECTID');

%%
%%calculate area of each polygon in raster%%
%is saved in polygon structure

%get unitarea
deg_km = 111.12;
%deg_km = 111.132;
dlat = abs(min(diff(lat(:,1)))); 
dlon = abs(min(diff(lon(1,:))));
dy = dlat*(deg_km * 1000);
dx = dlon*cosd(mean(lat(:,1)))*(deg_km * 1000); 
dx = 2.03326435575496;
dy = 2.97652555692817;
unitarea = dx*dy;
%dy dx are not exactly right and are different to the one when taking dx
%and dy in ArcGIS --> ArcGIS determines dx dy more precisely
%maybe use deg_km = 111.132;

%area for each slide
for i = 1:numel(S)
    raster_area = sum(sum(polyraster == S(i).OBJECTID))*unitarea;
    S(i).Rasterarea = raster_area;
    %difference Polygon (calculated in arcgis) and raster
    S(i).Areadiff = S(i).Shape_Area-S(i).Rasterarea;
end

%%
%%get starting and end points
%starting point: maximum of
start_area = zeros(size(polyraster));
end_area = zeros(size(polyraster));
length = zeros(size(polyraster));
for i = 1:numel(S)
    %extract slide
    %slide = polyraster == S(i).OBJECTID;
    slide = find(polyraster==S(i).OBJECTID);
    
    [~,imax] = max(elevation(slide));
    start_area(slide(imax)) = S(i).OBJECTID;
    
    %find maximum elevation of slide --> source area
    %[ymax,xmax] = find(elevation.*slide == max(elevation(slide)));
    %start_area(ymax,xmax) = S(i).OBJECTID;
    
    %find minimum elvation of slide --> end area
    [~,imin] = min(elevation(slide));
    end_area(slide(imin)) = S(i).OBJECTID;
    %[ymin,xmin] = find(elevation.*slide == min(elevation(slide)));
    %end_area(ymin,xmin) = S(i).OBJECTID;
    
    %get distance between start and end (with dx and dy)
    [ymax,xmax] = ind2sub(size(start_area),slide(imax));
    [ymin,xmin] = ind2sub(size(start_area),slide(imin));
    
    dys = abs(ymax-ymin)*dy;
    dxs = abs(xmax-xmin)*dx;
    dz = max(elevation(slide))-min(elevation(slide));
    lgt = sqrt(sqrt(dys^2+dxs^2)^2+dz^2);
    S(i).length = lgt;
end


%%
%%%%%ToDo%%%%
%*fläche berechnen --> zähle zellen für jeden slide(Ereignisnummer)
% vergleiche mit ARCGIS
%*source area finden --> höchster punkt auch mit alti3d2m (gleiche lat lon
% auflösung--> nicht gleiche x y auflösung)
%*flow distance berechenen --> distance höchster tiefster punkt (auch mit
%alti3d_2m)
%*bei berechnung von area wird slope nicht berücksichtigt
%*bei distanz wird slope nicht berücksichtig
%%%%%%%%%%%%%

end