function climada_ls_FPcalibration(calculate,calculate2)

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
%     calculate: 1 if area lenght start end need to be calcualted
%     calculate2: 1 if 
% OPTIONAL INPUT PARAMETERS:
%    
% OUTPUTS:
%    
% MODIFICATION HISTORY:
% Thomas Rölli, thomasroelli@gmail.com, 20180417, init

global climada_global
if ~climada_init_vars, return; end

if ~exist('calculate') calculate=0; end
if ~exist('calculate2') calculate2=0; end


load('C:\Users\Simon Rölli\Desktop\data\centroids_large\_LS_Sarnen_alti3d_2m.mat','centroids');
load('C:\Users\Simon Rölli\Desktop\data\centroids_large\_LS_Sarnen_interpol_alti3d.mat','centroids');

%read in DEM (original alti3d with 2x2m resolution is best choice --> best results when assessing source
%area lenght and width of slides
n_lon = numel(unique(centroids.lon));
n_lat = numel(unique(centroids.lat));
lon = reshape(centroids.lon,n_lat,n_lon);
lat = reshape(centroids.lat,n_lat,n_lon);
elevation = reshape(centroids.elevation_m,n_lat,n_lon);

n_lon_10m = numel(unique(centroids.lon));
n_lat_10m = numel(unique(centroids.lat));
lon_10m = reshape(centroids.lon,n_lat_10m,n_lon_10m);
lat_10m = reshape(centroids.lat,n_lat_10m,n_lon_10m);
elevation_10m = reshape(centroids.elevation_m,n_lat_10m,n_lon_10m);
clear centroids;

if calculate
    %read in shape file of past shallow landslides
    %S = shaperead('C:\Users\Simon Rölli\Desktop\data\inventory\slides_forMatlab.shp');
    S = shaperead('C:\Users\Simon Rölli\Desktop\data\inventory\be_polygon_forMatlab.shp');
    field = 'ID';

    [S,polyraster,start_area,end_area] = climada_ls_scoresShapefile(lon,lat,elevation,S,field,0);

%     figure
%     plot([S.AREA_GIS],[S.P_AREA_noSL],'*')
%     figure
%     plot([S.AREA_GIS],[S.R_AREA_noSL],'*')
%     figure
%     plot([S.AREA_GIS],[S.R_AREA_SL],'*')
%     figure
%     plot([S.AREA_GIS],[S.R_AREA_SL],'*')

    %calculate ratio of length/shape length
    C = num2cell([S.R_LGT_noSL]./[S.SHAPE_Leng]);
    [S.LGT_RAT_noSL] = C{:};
    %with slope
    C = num2cell([S.R_LGT_SL]./[S.SHAPE_Leng]);
    [S.LGT_RAT_SL] = C{:};

    shapewrite(S,'C:\Users\Simon Rölli\Desktop\data\calibration\be_fromMatlab.shp');
    grid2tiff(lon,lat,start_area,'C:\Users\Simon Rölli\Desktop\data\calibration\start_area.tif');
    grid2tiff(lon,lat,end_area,'C:\Users\Simon Rölli\Desktop\data\calibration\end_area.tif');
    grid2tiff(lon,lat,polyraster,'C:\Users\Simon Rölli\Desktop\data\calibration\polyraster.tif');
end
if calculate2
   S = shaperead('C:\Users\Simon Rölli\Desktop\data\calibration\be_fromMatlab.shp');
%    grid = GRIDobj('C:\Users\Simon Rölli\Desktop\data\calibration\polyraster.tif');
%    polyraster = flipud(grid.Z);
   grid = GRIDobj('C:\Users\Simon Rölli\Desktop\data\calibration\start_area.tif');
   start_area = flipud(grid.Z);
   grid = GRIDobj('C:\Users\Simon Rölli\Desktop\data\calibration\end_area.tif');
   end_area = flipud(grid.Z);
   clear grid;
   
    %load 10m interpolated alti3d data
    
    %find nearest cell in 10m grid and assign corresponding value 
    start_10m = zeros(n_lat_10m,n_lon_10m);
    end_10m = zeros(n_lat_10m,n_lon_10m);
    %prepare grid to use in dsearchn()
    grid_10m = [lon_10m(:) lat_10m(:)];
    
    %indices of start areas
    idx = find(start_area~=0);
    points = [lon(idx) lat(idx)];
    k = dsearchn(grid_10m,points);
    start_10m(k) = start_area(idx);
    
    %indices of end areas
    idx = find(end_area~=0);
    points = [lon(idx) lat(idx)];
    k = dsearchn(grid_10m,points);
    end_10m(k) = end_area(idx);
    
    grid2tiff(lon_10m,lat_10m,start_10m,'C:\Users\Simon Rölli\Desktop\data\calibration\start_10m.tif');
    grid2tiff(lon_10m,lat_10m,end_10m,'C:\Users\Simon Rölli\Desktop\data\calibration\end_10m.tif');
end

S = shaperead('C:\Users\Simon Rölli\Desktop\data\calibration\be_fromMatlab.shp');
grid = GRIDobj('C:\Users\Simon Rölli\Desktop\data\calibration\start_10m.tif');
start_10m = flipud(grid.Z);
grid = GRIDobj('C:\Users\Simon Rölli\Desktop\data\calibration\end_10m.tif');
end_10m = flipud(grid.Z);
clear grid;

dH = 0;
exponent = 25;
v_max = 4;
phi = 22;
delta_i = 0.0003;
perWt = [1 0.8 0.4 0 0 0 0.4 0.8];

source = zeros(size(elevation_10m));
elevation_10m = fillsinks(elevation_10m);
mult_flow = climada_ls_multipleflow(lon_10m,lat_10m,elevation_10m,exponent,dH,1);
[~,horDist,verDist] = climada_centroids_gradients(lon_10m,lat_10m,elevation_10m);

for i = 1:numel(S)
    source(find(start_10m==S(i).ID)) = 1;
    spreaded = climada_ls_propagation(source,mult_flow,horDist,verDist,v_max,phi,delta_i,perWt);
    
    
    source = source*0;
end


end