function raster = climada_ls_poly2raster(lon,lat,S,field,buffer)

% Transform polygons into raster data
% MODULE:
%   flood
% NAME:
%   climada_ls_poly2raster
% PURPOSE:
%   Transform polygons from a shape-structure into a raster-grid. Function
%   uses polygon2GRIDobj from TopoToolbox and adds an optional buffer. The
%   TopoToolbox function transforms each gridpoint enclosed by the polygon
%   (or lying within the polygon) to a rastercell and labels it according
%   to the name or number of the field of the polygon. When a buffer is set
%   the area is widened and accordinlgy, the number of cells lying within 
%   the polygon is increased.
% CALLING SEQUENCE:
%   climada_ls_hazard_sets
% EXAMPLE:
%   transformation with a buffer of 1arcsec. no field given --> logical
%   transformation (0/1).
%    raster = climada_ls_poly2raster(lon,lat,S,'',1/3600);
%    tranformation while it takes writes the information from TARGED_FID
%   for cells within the corresponding polygon. Without buffer
%    raster = climada_ls_poly2raster(lon,lat,S,'TARGET_FID');
% INPUTS: 
%    lon: (nxm)-matrix. longitutional information about the position of the desired grid
%    lat: (nxm)-matrix. latitutional information about the position of the desired grid
%    S:   mapstruct of a polyline as imported by shaperead. Must have
%         the fields X and Y with the coordinates of the 
% OPTIONAL INPUT PARAMETERS:
%    field: name of the numeric data to be mapped. If not
%           provided, polygon2GRIDobj will map logical values.
%    buffer: Scalar which is specified in degrees of
%            arc and defines the width of the bufferzone along the polygon
%            (default = 0)
% OUTPUTS:
%    raster: Transformed raster grid with same dimension as lat and lon
% MODIFICATION HISTORY:
% Thomas Rölli, thomasroelli@gmail.com, 20180327, init

global climada_global
if ~climada_init_vars, return; end

% check arguments
if ~exist('lon', 'var'), lon = []; end
if ~exist('lat', 'var'), lat = []; end
if ~exist('field', 'var'), field = []; end
if ~exist('buffer', 'var'), buffer = []; end

% PARAMETERS 
if isempty(lon); return; end
if isempty(lat); return; end
if isempty(field); field = []; end
if isempty(buffer); buffer = 0; end

DEM = GRIDobj(lon,lat,zeros(size(lon)));
if isempty(field); P = polygon2GRIDobj(DEM,S); else P = polygon2GRIDobj(DEM,S,field); end

%generate buffer around original shape file
if buffer > 0
    Sb = S;
    for i = 1:numel(S)
        [latb,lonb] = bufferm(S(i).Y,S(i).X,buffer);
        Sb(i).X = lonb;
        Sb(i).Y = latb;
    end
else
    Sb = [];
end


if isempty(field); Pb = polygon2GRIDobj(DEM,Sb); else Pb = polygon2GRIDobj(DEM,Sb,field); end

raster = P.Z+Pb.Z;
raster = flipud(raster);

end