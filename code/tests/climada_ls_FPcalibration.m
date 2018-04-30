function climada_ls_FPcalibration(calculate)

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

if calculate
    %read in shape file of past shallow landslides
    %S = shaperead('C:\Users\Simon Rölli\Desktop\data\inventory\slides_forMatlab.shp');
    S = shaperead('C:\Users\Simon Rölli\Desktop\data\inventory\be_polygon_forMatlab.shp');
    field = 'ID';

    %read in DEM (original alti3d with 2x2m resolution is best choice --> best results when assessing source
    %area lenght and width of slides
    load('C:\Users\Simon Rölli\Desktop\data\centroids_large\_LS_Sarnen_alti3d_2m.mat','centroids');

    n_lon = numel(unique(centroids.lon));
    n_lat = numel(unique(centroids.lat));
    lon = reshape(centroids.lon,n_lat,n_lon);
    lat = reshape(centroids.lat,n_lat,n_lon);
    elevation = reshape(centroids.elevation_m,n_lat,n_lon);

    [S,start_area,end_area] = climada_ls_scoresShapefile(lon,lat,elevation,S,field,0);

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
end

end