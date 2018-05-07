function [S,polyraster,src_ar,end_ar] = climada_calibration_orgS2raster(lon,lat,elevation,S,field,buffer)

% 
% MODULE:
%   flood
% NAME:
%   climada_ls_FPcalibration
% PURPOSE:
%   function transforms the shapefile with polygons into a raster with a
%   given resolution and derives scores of each individual polygon/slide
%   such as the area, length, source and end raster point.
% CALLING SEQUENCE:
%   
% EXAMPLE:
%   
% INPUTS: 
%   lon:       (nxm)-matrix. longitutional information about the position 
%              of the desired grid
%   lat:       (nxm)-matrix. latitutional information about the position 
%              of the desired grid
%   elevation: (nxm)-matrix. DEM corresponding to lon/lat matrices.
%              Resolution of DEM defines resolution of landslide raster
%              when transformed from polygons to raster.
%   S:         mapstruct of a polyline as imported by shaperead. Must have
%              the fields X and Y with the coordinates of the polygons (lying
%              withing the DEM).
%   field:     name (string, e.g. 'ID') of the numeric data to be mapped. Values in field need
%              to be numeric.
% OPTIONAL INPUT PARAMETERS:
%   buffer:    Scalar which is specified in degrees of
%              arc and defines the width of the bufferzone along the polygon
%              (default = 0) 
% OUTPUTS:
%   S:        mapstruct of a polyline as imported by shaperead. with additional fields:
%             .R_AREA_noSL: area of transformed raster of each polygon when
%              not considering the slope.
%             .R_AREA_SL: area of transformed raster of each polygon when
%              considering the slope.
%             .P_AREA_noSL: area of polygon using polyarea-function when
%              not considering the slope.
%             .R_LGT_noSL: length of transformed raster of each polygon
%              when not considering the slope. the lenght is defined as the
%              distance from the highest to the lowest raster point of the
%              corresponding polygon.
%             .R_LGT_noSL: length of transformed raster of each polygon
%              when considering the slope. the lenght is defined as the
%              distance from the highest to the lowest raster point of the
%              corresponding polygon.
%   polyraster: Transformed polygons in gridded raster labelled with the
%               corresponding name of the choosen 'field'
%   src_ar:   (nxm)-matrix with source area (highest points) of each slide.
%             The slides are labelled with the corresponding number given
%             in S.(field).
%   end_ar:   (nxm)-matrix with end area (lowest points) of each slide.
%             The slides are labelled with the corresponding number given
%             in S.(field).
% MODIFICATION HISTORY:
% Thomas Rölli, thomasroelli@gmail.com, 20180426, init

global climada_global
if ~climada_init_vars, return; end

% check arguments
if ~exist('lon', 'var'), lon = []; end
if ~exist('lat', 'var'), lat = []; end
if ~exist('elevation', 'var'), elevation = []; end
if ~exist('S', 'var'), S = []; end
if ~exist('field', 'var'), field = []; end
if ~exist('buffer', 'var'), buffer = []; end

% PARAMETERS 
if isempty(lon); return; end
if isempty(lat); return; end
if isempty(elevation); return; end
if isempty(S); return; end
if isempty(field); field = []; end
if isempty(buffer); buffer = 0; end

%init output
src_ar = zeros(size(elevation));
end_ar = zeros(size(elevation));

%%
%transform polygon in raster
polyraster = climada_ls_poly2raster(lon,lat,S,field);

%%
%%calculate area and length of each polygon in raster%%
%is saved in shape structure

%area of each raster cell
cell_area = climada_centroids_area(lon,lat,elevation,0);

%get unitarea
deg_km = 111.32;
%deg_km = 111.344; %best results are gotten when using polyarea()
dlat = abs(min(diff(lat(:,1)))); 
dlon = abs(min(diff(lon(1,:))));
dy = dlat*(deg_km * 1000);
dx = dlon*cosd(mean(lat(:,1)))*(deg_km * 1000); 
unitarea = dx*dy;

msgstr   = sprintf('Assigning properties to %i polygons ... ',numel(S));
mod_step = 10; % first time estimate after 10 assets, then every 100
if climada_global.waitbar
    fprintf('%s (updating waitbar with estimation of time remaining every 100th centroid)\n',msgstr);
    h        = waitbar(0,msgstr);
    set(h,'Name','Assigning TWI');
else
    fprintf('%s (waitbar suppressed)\n',msgstr);
    format_str='%s';
end

%iteration through all polygons --> derive area start end of polyraster
for i = 1:numel(S)
    %%
    %%%%area%%%%
    slide = find(polyraster==S(i).(field));
    
    %no consideration of slope --> with unitarea
    raster_area = numel(slide)*unitarea;
    %raster_area = sum(cell_area(slide));
    S(i).R_AREA_noSL = raster_area;
    
    %considers slope when calculating cell area
    raster_area_slope = sum(cell_area(slide));
    S(i).R_AREA_SL = raster_area_slope;
    
    %use polygon function of matlab to get area, therefore translate
    %lon lat to x y coordinates
    x = S(i).X*cosd(mean(lat(:,1)))*(deg_km * 1000);
    y = S(i).Y*deg_km*1000;
    %consider latitude dependency
    %meanLat = nanmean(S(i).Y);
    %x = S(i).X*cosd(meanLat)*(deg_km * 1000); %considers change with latitude
    poly_area = polyarea(x(1:numel(x)-1),y(1:numel(y)-1));
    S(i).P_AREA_noSL = poly_area;
    
    %%
    %%get starting (highest) and end (lowest) points
    %starting point: maximum of slide; end point: minimum of slide
    [~,imax] = max(elevation(slide));
    if ~isempty([lon(imax) lon(imax)]);
        S(i).startlon = lon(imax);
        S(i).startlat = lat(imax);
    else
        S(i).startlon = 0;
        S(i).startlat = 0;
    end
    src_ar(slide(imax)) = S(i).(field);
    
    %find minimum elvation of slide --> end area
    [~,imin] = min(elevation(slide));
    if ~isempty([lon(imin) lon(imin)])
        S(i).endlon = lon(imin);
        S(i).endlat = lat(imin);
    else
        S(i).endlon = 0;
        S(i).endlat = 0;
    end
    end_ar(slide(imin)) = S(i).(field);
    
    %get distance between start and end (with dx and dy)
    [ymax,xmax] = ind2sub(size(src_ar),slide(imax));
    [ymin,xmin] = ind2sub(size(src_ar),slide(imin));
    
    dys = abs(ymax-ymin)*dy;
    dxs = abs(xmax-xmin)*dx;
    dz = max(elevation(slide))-min(elevation(slide));
    try 
        lgt_noslope = double(sqrt(dys^2+dxs^2));
    catch
        lgt_noslope = double(0);
    end
    try 
        lgt = double(sqrt(sqrt(dys^2+dxs^2)^2+dz^2));
    catch
        lgt = double(0);
    end
    S(i).R_LGT_noSL = lgt_noslope;
    S(i).R_LGT_SL = lgt;
    
    if mod(i,mod_step)==0
        mod_step = 10;
        msgstr = sprintf('%i/%i polygons',i,numel(S));
        if climada_global.waitbar
            waitbar(i/numel(S),h,msgstr); % update waitbar
        else
            fprintf(format_str,msgstr); % write progress to stdout
            format_str=[repmat('\b',1,length(msgstr)) '%s']; % back to begin of line
        end
    end
end



end