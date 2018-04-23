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

%area of each slide
cell_area = climada_centroids_area(lon,lat,elevation);

%get unitarea
deg_km = 111.32;
deg_km = 111.344; %best results are gotten when using polyarea()
dlat = abs(min(diff(lat(:,1)))); 
dlon = abs(min(diff(lon(1,:))));
dy = dlat*(deg_km * 1000);
dx = dlon*cosd(mean(lat(:,1)))*(deg_km * 1000); 
dx_gis = 2.03326435575496;
dy_gis = 2.97652555692817;
unitarea = dx*dy;
unitarea_gis = dx_gis*dy_gis;
%dy dx are not exactly right and are different to the one when taking dx
%and dy in ArcGIS --> ArcGIS determines dx dy more precisely
%maybe use deg_km = 111.34;

for i = 1:numel(S)
    %%
    %%%%area%%%%
    slide = find(polyraster==S(i).OBJECTID);
    
    %no consideration of slope
    %matlab unitarea
    raster_area = numel(slide)*unitarea;
    S(i).matlab_unitarea = raster_area;
    %arcGis unitarea
    raster_area = numel(slide)*unitarea_gis;
    S(i).gis_unitarea = raster_area;
    
    %considers slope when calculating cell area
    %raster_area_slope = sum(cell_area(slide));
    %S(i).Matlab_area_slope = raster_area_slope;
    
    %difference Polygon (calculated in arcgis) and raster
    S(i).diff_mat = S(i).matlab_unitarea-S(i).Shape_Area;
    S(i).diff_gis = S(i).gis_unitarea-S(i).Shape_Area;
    %%
    %%get starting and end points
    %starting point: maximum of slide; end point: minimum of slide
    x = S(i).X*cosd(mean(lat(:,1)))*(deg_km * 1000);
    S(i).meanLat = nanmean(S(i).Y);
    
    %x = S(i).X*cosd(S(i).meanLat)*(deg_km * 1000);
    y = S(i).Y*deg_km*1000;
    
    poly_area = polyarea(x(1:numel(x)-1),y(1:numel(y)-1));
    S(i).matlab_polyarea = poly_area;
    S(i).diff_poly = S(i).matlab_polyarea-S(i).Shape_Area;
    
    
    
    i
end
% figure
% plot([S.Shape_Area],[S.diff_poly],'*')
% figure
% plot([S.meanLat],[S.diff_poly],'*')

figure
plot([S.Shape_Area],[S.diff_mat],'*')
figure
plot([S.Shape_Area],[S.diff_gis],'*')
figure
plot([S.Shape_Area],[S.diff_poly],'*')
figure
plot([S.Shape_Area],[S.diff_poly],'*')


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