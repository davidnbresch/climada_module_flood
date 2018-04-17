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

polyraster = climada_ls_poly2raster(lon,lat,S,'Ereignisnu');

[8.111086 8.341411 46.815369 46.946240]

end