function [IDfield,area,lgt,ID] = climada_ls_slideScores(intensity,source_ID,elevation,cell_area,dx,dy)

% Returns different scores of slides (length, area)
% MODULE:
%   flood
% NAME:
%   climada_ls_slideScores
% PURPOSE:
%   By providing start area with the source ID (source_ID) and the
%   intensity map after spreading (intensity), the function detects
%   connected area and labels the the areas with the corresponding ID.
%   For each individual slide the area and length is calculated.
% CALLING SEQUENCE:
%   climada_ls_slideScores(intensity,source_ID,elevation,cell_area,dx,dy)
% EXAMPLE:
%   [IDfield,area,lgt,ID] = climada_ls_slideScores(intensity,source_ID,elevation,cell_area,dx,dy)
% INPUTS: 
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
%    area:      1xs vector with area of each individual slide (s --> number of slides in source_ID).
%    lgt:       1xs vector with length of each individual slide. Length of
%               straight line from start to end (lowest point) of slide by
%               assuming constant linear slope.
%    ID:        1xs vector with slide id of each corresponding slide given in source_ID.
% MODIFICATION HISTORY:
% Thomas Rölli, thomasroelli@gmail.com, 20180502, init

IDfield = zeros(size(intensity));

%get ID number of corresponding spreaded slide --> connect areas and assign
%ID of source cell to the whole area.
con_areas = bwconncomp(intensity>0); %detection of connected areas
con_areas_idx = con_areas.PixelIdxList;
source_idx = find(source_ID~=0)';
%find fid in cell structure --> gives out position
c = cellfun(@(x)(ismember(source_idx,x)),con_areas_idx,'UniformOutput',false);
[aidx,~] = find(reshape([c{:}],numel(source_idx),[])'); %find position of 1 in c --> index gives corresponding position of indices in con_areas_idx

area = zeros(size(source_idx));
lgt = zeros(size(source_idx));
ID = zeros(size(source_idx));

for n=1:numel(source_idx)
    slide_idx = cell2mat(con_areas_idx(aidx(n)));
    slide_id = source_ID(source_idx(n));
    ID(n) = slide_id;
    IDfield(slide_idx) = slide_id;
    area(n) = sum(cell_area(slide_idx));
    
    %calculate distance of slide in new raster (with slope)
    [Ist,Jst] = ind2sub(size(intensity),source_idx(n));
    [~,end_idx] = min(elevation(slide_idx));
    [Ien,Jen] = ind2sub(size(intensity),slide_idx(end_idx));
    dI = abs(Ist-Ien)*dy;
    dJ = abs(Jst-Jen)*dx;
    dz = elevation(source_idx(n))-elevation(slide_idx(end_idx));
    try 
        lgt(n) = double(sqrt(sqrt(dI^2+dJ^2)^2+dz^2));
    catch
        lgt(n) = double(0);
    end   
end

end