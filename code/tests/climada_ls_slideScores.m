function [IDfield,S] = climada_ls_slideScores(lon,lat,intensity,source_ID,elevation,cell_area,field)

% Returns different scores of slides (length, area)
% MODULE:
%   flood
% NAME:
%   climada_ls_slideScores
% PURPOSE:
%   By providing start area with the source ID (source_ID) and the
%   intensity map after spreading (intensity), the function detects
%   connected area and labels the the areas with the corresponding ID.
%   For each individual slide the area and length is calculated. It is
%   important that the slides are located such that slides do not flow over
%   the same areas --> will result in too large areas
% CALLING SEQUENCE:
%   climada_ls_slideScores(intensity,source_ID,elevation,cell_area,dx,dy)
% EXAMPLE:
%   [IDfield,area,lgt,ID] = climada_ls_slideScores(intensity,source_ID,elevation,cell_area,dx,dy)
% INPUTS:
%     lon/lat:   nxm matrix with longitudinal/latitudinal information of
%                grid
%     intensity: nxm matrix with the spreaded intensity of each considered
%                slide. The affected areas need to be seperated by each
%                other.
%     source_ID: nxm matrix with corresponding source area ID. The source
%                areas need to correspond to the slides in intensity
%     elevation: nxm matrix with elevation data --> needed to find end area
%               (lowest point)
%     cell_area: nxm matrix with area of each cell. calculated by:
%                climada_centroids_area
%     dx:        longitudinal resolution in [m], assumed constant within
%                study area.
%     dy:        latitudinal resolution in [m] 
% OPTIONAL INPUT PARAMETERS:
%    
% OUTPUTS:
%    IDfield:   nxm matrix with same areas as intensity but labelled with
%               corresponding id
%    S          Structure with:
%               .area: 1xs vector with area of each individual slide (s --> number of slides in source_ID).
%               .lgt: 1xs vector with length of each individual slide. Length of
%               straight line from start to end (lowest point) of slide by
%               assuming constant linear slope.
%               .ID: 1xs vector with slide id of each corresponding slide given in source_ID.
%               .dL: horizontal distance between start and end -> slope =0
%               .dH: vertical distance between start and end
%               
% MODIFICATION HISTORY:
% Thomas Rölli, thomasroelli@gmail.com, 20180502, init

%check arguments
if ~exist('lon'), return; end
if ~exist('lat'), return; end
if ~exist('intensity'), return; end
if ~exist('source_ID'), return; end
if ~exist('elevation'), return; end
if ~exist('cell_area'), return; end
if ~exist('field'), field = 'ID'; end

%init output
IDfield = zeros(size(intensity));
S = [];

%calculate dx and dy
deg_km = 111.32;
dlat = abs(min(diff(lat(:,1)))); 
dlon = abs(min(diff(lon(1,:))));
dy = dlat*(deg_km * 1000);
dx = dlon*cosd(mean(lat(:,1)))*(deg_km * 1000);

%get ID number of corresponding spreaded slide --> connect areas and assign
%ID of source cell to the whole area.
con_areas = bwconncomp(intensity>0); %detection of connected areas
con_areas_idx = con_areas.PixelIdxList;
source_idx = find(source_ID~=0)';
%find fid in cell structure --> gives out position
c = cellfun(@(x)(ismember(source_idx,x)),con_areas_idx,'UniformOutput',false);
[aidx,~] = find(reshape([c{:}],numel(source_idx),[])'); %find position of 1 in c --> index gives corresponding position of indices in con_areas_idx

%init fields which will be saved in S.field
startlon = zeros(size(source_idx));
endlon = zeros(size(source_idx));
startlat = zeros(size(source_idx));
endlat = zeros(size(source_idx));
area = zeros(size(source_idx));
lgt = zeros(size(source_idx));
ID = zeros(size(source_idx));
dH = zeros(size(source_idx));
dL = zeros(size(source_idx));

for n=1:numel(source_idx)
    slide_idx = cell2mat(con_areas_idx(aidx(n)));
    
    %write ID
    slide_id = source_ID(source_idx(n));
    
    ID(n) = slide_id;

    IDfield(slide_idx) = slide_id;
    
    %write area
    area(n) = sum(cell_area(slide_idx));
    
    %calculate distance of slide in new raster (with slope)
    [Ist,Jst] = ind2sub(size(intensity),source_idx(n));
    [~,end_idx] = min(elevation(slide_idx));
    [Ien,Jen] = ind2sub(size(intensity),slide_idx(end_idx));
    dI = abs(Ist-Ien)*dy;
    dJ = abs(Jst-Jen)*dx;
    dz = elevation(source_idx(n))-elevation(slide_idx(end_idx));
    try 
        dL(n) = double(sqrt(dI^2+dJ^2));
        lgt(n) = double(sqrt(dL^2+dz^2));
    catch
        dL(n) = double(0);
        lgt(n) = double(0);
    end
    
    %write lon/lat of start and end points
    dH(n) = dz;
    startlon(n) = lon(source_idx(n));
    startlat(n) = lat(source_idx(n));
    endlon(n) = lon(slide_idx(end_idx));
    endlat(n) = lat(slide_idx(end_idx));
    
end

for i = 1:numel(startlon)
   S(i).startlon = startlon(i); 
   S(i).endlon = endlon(i); 
   S(i).startlat = startlat(i); 
   S(i).endlat = endlat(i); 
   S(i).(field) = ID(i); 
   S(i).area = area(i); 
   S(i).length = lgt(i); 
   S(i).dL = dL(i);
   S(i).dH = dH(i);
end



end