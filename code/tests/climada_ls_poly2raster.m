function raster = climada_ls_poly2raster(lon,lat,S,field,buffer)

% Transform polygons into raster data
% MODULE:
%   flood
% NAME:
%   climada_ls_poly2raster
% PURPOSE:
%   Transform polygons from a shape-structure into a raster-grid. Function
%   uses polygon2GRIDobj from TopoToolbox and adds an optional buffer. The
%   buffer is useful when the resolution of the DEM is sparse (without
%   buffer: would be inconsistent with too small area.)
% CALLING SEQUENCE:
%   climada_ls_hazard_sets
% EXAMPLE:
%   
% INPUTS: 
%    
% OPTIONAL INPUT PARAMETERS:
% 
% OUTPUTS:
%  
% MODIFICATION HISTORY:
% Thomas Rölli, thomasroelli@gmail.com, 20180327, init

global climada_global
if ~climada_init_vars, return; end

% check arguments
if ~exist('lon', 'var'), lon = []; end
if ~exist('lat', 'var'), lat = []; end
if ~exist('polygon', 'var'), polygon = []; end
if ~exist('field', 'var'), end
if ~exist('buffer', 'var'), buffer = []; end

% PARAMETERS 
if isempty(lon); return; end
if isempty(lat); return; end
if isempty(polygon); return; end
if isempty(field); end
if isempty(buffer); buffer = 0; end

DEM = GRIDobj(lon,lat,zeros(size(lon)));
P = polygon2GRIDobj(DEM,S,'TARGET_FID');

end