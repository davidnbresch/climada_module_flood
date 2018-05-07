function climada_ls_FPcalibration(calculate,calculate2)

% Script to calibrate flow path parameters
% INPUTS: 
%     calculate: 1 if area lenght start end of original high resolution Sample need to be calcualted
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


%read in DEM (original alti3d with 2x3m resolution is best choice --> best results when assessing source
%area lenght and width of slides 
orgSave = '2x3m';
load('C:\Users\Simon Rölli\Desktop\data\centroids_large\_LS_Sarnen_alti3d_2m.mat','centroids');
n_lon = numel(unique(centroids.lon));
n_lat = numel(unique(centroids.lat));
orgLon = reshape(centroids.lon,n_lat,n_lon);
orgLat = reshape(centroids.lat,n_lat,n_lon);
orgElevation = reshape(centroids.elevation_m,n_lat,n_lon);

calSave = '7x10m';
%load('C:\Users\Simon Rölli\Desktop\data\centroids_large\_LS_Sarnen_srtm1','centroids');
load('C:\Users\Simon Rölli\Desktop\data\centroids_large\_LS_Sarnen_interpol_alti3d','centroids');
n_lon = numel(unique(centroids.lon));
n_lat = numel(unique(centroids.lat));
lon = reshape(centroids.lon,n_lat,n_lon);
lat = reshape(centroids.lat,n_lat,n_lon);
elevation = reshape(centroids.elevation_m,n_lat,n_lon);
clear centroids;

field = 'ID';

if calculate
    %read in shape file of past shallow landslides
    %S = shaperead('C:\Users\Simon Rölli\Desktop\data\inventory\slides_forMatlab.shp');
    orgS = shaperead('C:\Users\Simon Rölli\Desktop\data\inventory\be_polygon_forMatlab.shp');

    [orgS,orgPolyraster,orgStart,orgEnd] = climada_calibration_orgS2raster(orgLon,orgLat,orgElevation,orgS,field,0);

%     figure
%     plot([S.AREA_GIS],[S.P_AREA_noSL],'*')
%     figure
%     plot([S.AREA_GIS],[S.R_AREA_noSL],'*')
%     figure
%     plot([S.AREA_GIS],[S.R_AREA_SL],'*')
%     figure
%     plot([S.AREA_GIS],[S.R_AREA_SL],'*')
%     figure
%     plot([orgS.AREA_GIS],([orgS.R_AREA_noSL]-[orgS.AREA_GIS])./[orgS.AREA_GIS],'*')
%     xlabel('ArcGIS [m^2]')
%     ylabel('rel. Difference Matlab raster vs ArcGIS [%]')

    %calculate ratio of length/shape length
    C = num2cell([orgS.R_LGT_noSL]./[orgS.SHAPE_Leng]);
    [orgS.LGT_RAT_noSL] = C{:};
    %with slope
    C = num2cell([orgS.R_LGT_SL]./[orgS.SHAPE_Leng]);
    [orgS.LGT_RAT_SL] = C{:};

    shapewrite(orgS,['C:\Users\Simon Rölli\Desktop\data\calibration\shapes\orgS_' orgSave '.shp']);
    grid2tiff(orgLon,orgLat,orgStart,['C:\Users\Simon Rölli\Desktop\data\calibration\orgStart_' orgSave '.tif']);
    grid2tiff(orgLon,orgLat,orgEnd,['C:\Users\Simon Rölli\Desktop\data\calibration\orgEnd_' orgSave '.tif']);
    grid2tiff(orgLon,orgLat,orgPolyraster,['C:\Users\Simon Rölli\Desktop\data\calibration\orgPolyraster_' orgSave '.tif']);
end
if calculate2
   orgS = shaperead(['C:\Users\Simon Rölli\Desktop\data\calibration\shapes\orgS_' orgSave '.shp']);
   grid = GRIDobj(['C:\Users\Simon Rölli\Desktop\data\calibration\orgStart_' orgSave '.tif']);
   orgStart = flipud(grid.Z);
   grid = GRIDobj(['C:\Users\Simon Rölli\Desktop\data\calibration\orgEnd_' orgSave '.tif']);
   orgEnd = flipud(grid.Z);
   clear grid;
   
   %get length and starting/end points of slides when changing to coarser
   %resolution (snapping)
   [snapS,snapStart,snapEnd_] = climada_calibration_snap(orgS,orgLon,orgLat,orgStart,orgEnd,lon,lat,elevation,field);
   
   %%
   %%%toDo
   %get endlon startlon to X and lat to Y --> to be able to draw line in
   %arcgis
   
   %%
   
   shapewrite(orgS,['C:\Users\Simon Rölli\Desktop\data\calibration\shapes\snapS_' calSave '.shp']);
   grid2tiff(orgLon,orgLat,orgStart,['C:\Users\Simon Rölli\Desktop\data\calibration\snapStart_' calSave '.tif']);
   grid2tiff(orgLon,orgLat,orgEnd,['C:\Users\Simon Rölli\Desktop\data\calibration\snapEnd_' calSave '.tif']);
   
   
   [S,spreaded,IDfield] = climada_ls_FPcalibration(orgS,orgLon,orgLat,orgStart,orgEnd,lon,lat,elevation,field);
    
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
    
    grid2tiff(lon_10m,lat_10m,start_10m,'C:\Users\Simon Rölli\Desktop\data\calibration\start_30m.tif');
    grid2tiff(lon_10m,lat_10m,end_10m,'C:\Users\Simon Rölli\Desktop\data\calibration\end_30m.tif');
end




end