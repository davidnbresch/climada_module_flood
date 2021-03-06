function [subS,polyraster,src_ar,end_ar] = climada_calibration_orgS2raster(lon,lat,elevation,S,field,buffer)

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
%   subS:      mapstruct of a polgones as imported by shaperead. with additional fields:
%             .R_AREA_noSL: area of transformed raster of each polygon when
%              not considering the slope.
%             .area: area of transformed raster of each polygon when
%              considering the slope.
%             .P_AREA_noSL: area of polygon using polyarea-function when
%              not considering the slope.
%             .dL: length of transformed raster of each polygon
%              when not considering the slope (horizontal length of slide). the lenght is defined as the
%              horizontal distance from the highest to the lowest raster point of the
%              corresponding polygon.
%             .length: length of transformed raster of each polygon
%              when considering the slope. the lenght is defined as the
%              distance from the highest to the lowest raster point of the
%              corresponding polygon.
%             .dH: vertical difference of start (highest) and end cell
%             (lowest)
%             .removed slides which are not covered by the high resolution
%              grid. Most probably because they are too small.
%             .slope saves slope at starting cell (slope from
%             climada_centroids_slope)
%             .reachAngle saves slope of line from start to end cell,
%             corresponds to angle of reach
%   subS:      mapstruct of polylines with fields:
%             .ID
%             .X/.Y
%             .length
%             .area
%   polyraster: Transformed polygons in gridded raster labelled with the
%               corresponding name of the choosen 'field'
%   src_ar:   (nxm)-matrix with source area (highest points) of each slide.
%             The slides are labelled with the corresponding number given
%             in S.(field).
%   end_ar:   (nxm)-matrix with end area (lowest points) of each slide.
%             The slides are labelled with the corresponding number given
%             in S.(field).
% MODIFICATION HISTORY:
% Thomas R�lli, thomasroelli@gmail.com, 20180426, init
% Thomas R�lli, thomasroelli@gmail.com, 20180522, add start slope and
%  slide slope

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

%slope of each raster cell
slope = climada_centroids_slope(lon,lat,elevation);

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
    
    subS(i).(field) = S(i).(field);
    
    %no consideration of slope --> with unitarea
    raster_area = numel(slide)*unitarea;
    %raster_area = sum(cell_area(slide));
    subS(i).R_AREA_noSL = raster_area;

    %mark slides which are not captured by the grid
    if(raster_area == 0)
        subS(i).removed = 1;
    else
        subS(i).removed = 0;
    end
    
    %considers slope when calculating cell area
    raster_area_slope = sum(cell_area(slide));
    subS(i).area = raster_area_slope;
    
    %use polygon function of matlab to get area, therefore translate
    %lon lat to x y coordinates
    x = S(i).X*cosd(mean(lat(:,1)))*(deg_km * 1000);
    y = S(i).Y*deg_km*1000;
    %consider latitude dependency
    %meanLat = nanmean(S(i).Y);
    %x = S(i).X*cosd(meanLat)*(deg_km * 1000); %considers change with latitude
    poly_area = polyarea(x(1:numel(x)-1),y(1:numel(y)-1));
    subS(i).P_AREA_noSL = poly_area;
    
    %%
    %%get starting (highest) and end (lowest) points
    %starting point: maximum of slide; end point: minimum of slide
    [~,imax] = max(elevation(slide));
    if ~isempty([lon(slide(imax)) lon(slide(imax))]);
        subS(i).startlon = lon(slide(imax));
        subS(i).startlat = lat(slide(imax));
    else
        subS(i).startlon = 0;
        subS(i).startlat = 0;
    end
    src_ar(slide(imax)) = S(i).(field);
    
    %find minimum elvation of slide --> end area
    [~,imin] = min(elevation(slide));
    if ~isempty([lon(slide(imin)) lon(slide(imin))])
        subS(i).endlon = lon(slide(imin));
        subS(i).endlat = lat(slide(imin));
    else
        subS(i).endlon = 0;
        subS(i).endlat = 0;
    end
    end_ar(slide(imin)) = S(i).(field);
    
    %get distance between start and end (with dx and dy)
    [ymax,xmax] = ind2sub(size(src_ar),slide(imax));
    [ymin,xmin] = ind2sub(size(src_ar),slide(imin));
    
    dys = abs(ymax-ymin)*dy;
    dxs = abs(xmax-xmin)*dx;
    dz = max(elevation(slide))-min(elevation(slide));
    %horizontal distance
    try 
        lgt_noslope = double(sqrt(dys^2+dxs^2));
    catch
        lgt_noslope = double(0);
    end
    %length with slope
    try 
        lgt = double(sqrt(sqrt(dys^2+dxs^2)^2+dz^2));
    catch
        lgt = double(0);
    end
    subS(i).dL = lgt_noslope;
    subS(i).length = lgt;
    subS(i).dH = double(dz);
    
    %write slope of start and slide (from start to end) in S.
    subS(i).slope = slope(ymax,xmax);
    try 
        subS(i).reachAngle = double(asind(dz/lgt));
    catch
        subS(i).reachAngle = double(0);
    end
    if(raster_area == 0)
        subS(i).slope = double(0);
        subS(i).reachAngle = double(0);
        subS(i).dH = double(0);
    end
    
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

%create new shape file with only most important info --> with polyline XY
%coordinates
%write polyline (X,Y) coordinates in structure, together with ID, area and length
X = [subS.startlon;subS.endlon]';
Y = [subS.startlat;subS.endlat]';

%set coordinates of removed slides to nan
r_idx = find([subS.removed]~=0);
X(r_idx,:) = nan;
Y(r_idx,:) = nan;
% dum = [snapS.(field)];
% ID(r_idx) = dum(r_idx);


for i = 1:numel(X(:,1))
    subS(i).Geometry = 'Line';
    subS(i).X = X(i,:);
    subS(i).Y = Y(i,:);
end

subS = rmfield(subS,{'startlon';'startlat';'endlon';'endlat'})

end