function [S,start,end_] = climada_calibration_snap(orgS,orgLon,orgLat,orgStart,orgEnd,...
    lon,lat,elevation,field)

% 
% MODULE:
%   flood
% NAME:
%   climada_calibration_snap
% PURPOSE:
%   Transforms starting and end points of slides of original dataset
%   (indicated in orgStart/orgEnd --> derived from orgS) to new grid with
%   coarser resolution of lon/lat grid. Therefore it uses the dsearchn
%   function to find the nearest grid point in lon/lat of each
%   starting/ending point in orgStart/orgEnd. If several starting points
%   are assigned to the same cell in the new/coarser grid, the starting
%   point is assigned to the slide with the greater area (in orgS). The
%   left out slides saved with []/0 in S.
% CALLING SEQUENCE:
%   [S,start,end_] = climada_calibration_snap(orgS,orgLon,orgLat,orgStart,orgEnd,...
%   lon,lat,field)
% EXAMPLE:
%   
% INPUTS: 
%   orgS:      original S --> mapstruct of a polyline as imported by shaperead
%              with additional fields calculated in climada_calibration_orgS2raster
%              Must at least include field definded in 'field' (e.g. 'ID'),
%              in order to be able to compare derived slide scores (from
%              this function) with original S.
%
%   orgLat/lat: norgxmorg Matrix of the latitutional information about the grid cells 
%               corresponding to the original/coarser resolved slides.
%   orgLon/lon: nxm Matrix of the longitudinal information about the grid cells 
%               corresponding to the original/coarser resolved slides.
%   orgStart:   original starting points of the slides indicated with
%               'field'. Must have same dimension as orgLat/orgLon(norgxmorg)
%   orgEnd:     original end points of the slides indicated with
%               'field'. Must have same dimension as orgLat/orgLon (norgxmorg)
%   elevation:  nxm matrix with elevation data --> needed to find end area
%               (lowest point)
%   field:     name (string, e.g. 'ID') of the numeric data to be mapped. Values in field need
%              to be numeric.
% OPTIONAL INPUT PARAMETERS:
%  
% OUTPUTS:
%   S:        mapstruct of a Line.
%             .X: X vertices of lines from start(highest) to end(lowest) point of slide.
%             structure is such that they can be saved and plotted as a shapefile
%             .Y: Y vertices...
%             .removed: info about slides which were removed (=2) or
%             replaced by longer slide (=1)
%             .(field): label of slide (field = 'ID' is recommended). .(field)
%             must be a field in orgS
%             .snap_length: length of slides after they are snapped to the
%             new coarser grid
%   start:   (nxm)-matrix with source area (highest points) of each slide.
%             The slides are labelled with the corresponding number given
%             in S.(field). If several slides with same start cell --> ID
%             of longer slide is saved in matrix.
%   end_:    (nxm)-matrix with end area (lowest points) of each slide.
%             The slides are labelled with the corresponding number given
%             in S.(field). If several slides with same end cell --> one
%             which is processed later is saved --> random, therefore not
%             recommended to compare end_ar of single slides
% MODIFICATION HISTORY:
% Thomas R�lli, thomasroelli@gmail.com, 20180426, init
% Thomas R�lli, thomasroelli@gmail.com, 20180514, treatment of replaced
%  slides

global climada_global
if ~climada_init_vars, return; end

%init S
S = [];

n_lon = size(lon,2);
n_lat = size(lat,1);
start = zeros(n_lat,n_lon);
end_ = zeros(n_lat,n_lon);
removed = [orgS.removed]; %saves ID if different slide is taken (when two slides start at same cell after snapping -> takes longer

%caculation of unitarea --> not considering slope
deg_km = 111.32;
dlat = abs(min(diff(lat(:,1)))); 
dlon = abs(min(diff(lon(1,:))));
dy = dlat*(deg_km * 1000);
dx = dlon*cosd(mean(lat(:,1)))*(deg_km * 1000); 
unitarea = dx*dy;

%prepare grid to use in dsearchn()
grid = [lon(:) lat(:)];

%when get to coarser resolution --> 
%indices of start areas
fprintf('Find start cells of new/coarser grid\n');
idx_st = find(orgStart~=0); %indices of original start
points = [orgLon(idx_st) orgLat(idx_st)];
k_st = dsearchn(grid,points); 

fprintf('Find end cells of new/coarser grid\n');
idx_en = find(orgEnd~=0); %indices of original end
points = [orgLon(idx_en) orgLat(idx_en)]; %original points
k_en = dsearchn(grid,points); %get indices of nearest point to original points

maxfield='AREA_GIS';
start = zeros(size(lon));
end_ = zeros(size(lon));
length = zeros(size(lon)); %to check if other (longer) slide exists at same cell already
for i=1:numel(orgS)
    if i==646
      disp('')  
    end
    %if slide was removed before then skip
    S(i).(field) = orgS(i).(field); %transform field ID
    if ~(orgS(i).removed)
%     if i==1080
%         disp('')
%     end
    ik_st = find(orgStart(idx_st) == orgS(i).(field)); %find corresponding index of field ID
    ik_en = find(orgEnd(idx_en) == orgS(i).(field));
    %write start and end coordinates in Structure
    S(i).startlon = lon(k_st(ik_st));
    S(i).startlat = lat(k_st(ik_st));
    S(i).endlon = lon(k_en(ik_en));
    S(i).endlat = lat(k_en(ik_en));
    
    if orgS(i).R_AREA_noSL == 0
        S(i).startlon = nan;
        S(i).startlat = nan;
        S(i).endlon = nan;
        S(i).endlat = nan;
    end
    
    %calculate distance of slide in new raster (with slope)
    [Ist,Jst] = ind2sub(size(lon),k_st(ik_st));
    [Ien,Jen] = ind2sub(size(lon),k_en(ik_en));
    dI = abs(Ist-Ien)*dy;
    dJ = abs(Jst-Jen)*dx;
    dz = elevation(k_st(ik_st))-elevation(k_en(ik_en));
    try 
        lgt = double(sqrt(sqrt(dI^2+dJ^2)^2+dz^2));
    catch
        lgt = double(0);
    end
    %write snapped length in structure
    S(i).snap_length = lgt;
    
    %create start matrix --> if one cell is affected several times-->
    %take the one which is longer (of original raster length)
    %for end matrix not so important --> just overwrite
    old_lgt = length(k_st(ik_st));
    if (old_lgt == 0) ||(old_lgt < orgS(i).R_LGT_SL)
        length(k_st(ik_st)) = orgS(i).R_LGT_SL;
        start(k_st(ik_st)) = orgS(i).ID;
        %removed(i) = 0;
    else 
        %removed(i) = 1;
    end
    
    end_(k_en(ik_en)) = orgS(i).ID;
    
    else
    %for removed slides
    S(i).startlon = nan;S(i).startlat = nan;
    S(i).endlon = nan;S(i).endlat = nan;
    S(i).snap_length = nan;
    end %end of skip removed
    
end 


%write X and Y coordinates in new S such that it can be plotted as polyline
X = [S.startlon;S.endlon]';
Y = [S.startlat;S.endlat]';

for i = 1:numel(S)
    k(i).Geometry = 'Line';
    k(i).X = X(i,:);
    k(i).Y = Y(i,:);
end

%find removed (2) or replaced (1) slides and save info in S.removed
r_idx = ~ismember([S.(field)],unique(start));
removed(r_idx) = removed(r_idx)+1;
m = num2cell(removed);
[k.removed] = m{:};

m = num2cell([S.(field)]);
[k.(field)] = m{:};
m = num2cell([S.snap_length]);
[k.length] = m{:};
m = num2cell([S.snap_length]);
[k.length] = m{:};

S = k;


end
