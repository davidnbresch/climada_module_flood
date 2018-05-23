function [S,spreaded,IDfield] = climada_calibration_propagate(lon,lat,elevation,flowPara,cell_area,snapStart,snapS,field)

% Script to calibrate flow path parameters
% MODULE:
%   flood
% NAME:
%   climada_ls_FPcalibration
% PURPOSE:
%   Given the starting points (snapStart from snapS) and set of flow path
%   parameters, the slides are triggered and propagated. The length and
%   area is calculated by climada_ls_slideScores and saved in S with the
%   corresponding ID.
% CALLING SEQUENCE:
%   
% EXAMPLE:
%   
% INPUTS: 
%     lon/lat:   nxm matrix with longitudinal/latitudinal information of
%                grid (same resolution as used for snapS)
%     elevation: (nxm)-matrix. DEM corresponding to lon/lat matrices.
%     flowPara:  structure of flow parameters which are used by the flow
%                path algorithm. Structure need to include the fields with
%                exactly the same names:
%                .dH --> elevation of center cell, used in climada_ls_multipleflow
%                .exp --> exponent of Holmgren algorithm, used in climada_ls_multipleflow
%                .phi --> minimum travel angle or Fahrböschung, used in
%                climada_ls_propagation
%                .vmax --> maximum velocity, used in climada_ls_propagation
%                .iT --> intensity threshold, used in climada_ls_propagation
%                .perWT --> persistence weights vector (8 elements), used in climada_ls_propagation
%     cell_area: nxm matrix with area of each cell. calculated by:
%                climada_centroids_area
%     snapStart: (nxm)-matrix with source area (highest points) of each slide.
%                output of climada_calibration_snap. Slides are triggered
%                from the points given in this matrix.
%     snapS:     mapstruct of a Line. given by: climada_calibration_snap
%     field:     name (string, recommended: 'ID') of the numeric data to be mapped. Values in field need
%                to be numeric.
% OPTIONAL INPUT PARAMETERS:
%    
% OUTPUTS:
%     S:        mapstruct of a Line of slides after triggered and propagated by flow
%               path algorithm
%               .X: X vertices of lines from start(highest) to end(lowest) point of slide.
%               Was defined by taking corresponding coordinates of start
%               grid in snapStart.
%               structure is such that they can be saved and plotted as a shapefile
%               .Y: Y vertices...
%               .removed: info about slides which were removed (=2) or
%               replaced by longer slide (=1)
%               .(field): label of slide (field = 'ID' is recommended). .(field)
%               must be a field in snapS
%               .length: length of modelled slides after they are triggered
%               in given DEM
%               .area: area of modelled slides after they were triggered
%     spreaded: (mxn) matrix of intensity of modelled slides after slides were triggered
%               propagated starting from fixed source areas. In regions
%               with several slides maximum intensity is taken.
%     IDfield:  (mxn) matrix of ID (defined in S.field) of slides --> cells
%               which is affected by e.g. slide with ID = 1 are labelled
%               as 1. If several slides affect the same are --> the minimum
%               of the ID is taken.
% MODIFICATION HISTORY:
% Thomas Rölli, thomasroelli@gmail.com, 20180502, init
% Thomas Rölli, thomasroelli@gmail.com, 20180514, processing managment

global climada_global
if ~climada_init_vars, return; end

%check arguments
if ~exist('lon'), return; end
if ~exist('lat'), return; end
if ~exist('elevation'), return; end
if ~exist('flowPara'), return; end
if ~exist('cell_area'), return; end
if ~exist('snapStart'), return; end
if ~exist('snapS'), return; end
if ~exist('field'), field='ID'; end

%check flow path parameters
if ~isfield(flowPara,'dH'), dH = 0; else, dH = flowPara.dH; end
if ~isfield(flowPara,'exp'), exp = 25; else, exp = flowPara.exp; end
if ~isfield(flowPara,'phi'), phi = 22; else, phi = flowPara.phi; end
if ~isfield(flowPara,'vmax'), vmax = 4; else, vmax = flowPara.vmax; end
if ~isfield(flowPara,'iT'), iT = 0.0003; else, iT = flowPara.iT; end
if ~isfield(flowPara,'perWT'), perWT = [1 0.8 0.4 0 0 0 0.4 0.8]; else, perWT = flowPara.perWT; end

%init output
S = [];

n_lon = size(lon,2);
n_lat = size(lat,1);

%caculation of unitarea --> not considering slope
deg_km = 111.32;
dlat = abs(min(diff(lat(:,1)))); 
dlon = abs(min(diff(lon(1,:))));
dy = dlat*(deg_km * 1000);
dx = dlon*cosd(mean(lat(:,1)))*(deg_km * 1000); 
unitarea = dx*dy;


%needed data for flowpath
source = zeros(size(elevation));
elevation = fillsinks(elevation);
mult_flow = climada_ls_multipleflow(lon,lat,elevation,exp,dH,1);
[~,horDist,verDist] = climada_centroids_gradients(lon,lat,elevation);

%init matrices to save results
source = zeros(size(elevation));
source(find(snapStart~=0)) = 1;
sel_source = zeros(size(elevation));
mask = zeros(size(elevation));
spreaded = zeros(size(elevation)); %saves intensity (max when several slides)
temp_spreaded = zeros(size(elevation));
IDfield = zeros(size(elevation))+100000; %saves ID of slides (min when several slides --> +100000

%init fields in S
startlon = zeros(size(S));
endlon = zeros(size(S));
startlat = zeros(size(S));
endlat = zeros(size(S));
ID = zeros(size(S));
lgt = zeros(size(S));
area = zeros(size(S));

%%
%method which chooses several slide sources at a time to save computing
%time --> is done by taking slides which are far enough apart from each
%other (distance defined by buf_m)
numsource = numel(find(source==1)); %for process management
format_str = '%s';

buf_m = 1000; %bufferregion in meters in which no other slides are choosen--> prevents slides from flow over each other
imask = ceil(buf_m/dy);
jmask = ceil(buf_m/dx);
while sum(source(:)) > 0
    for i=1:n_lat
        for j=1:n_lon
            if (source(i,j) == 1) && (mask(i,j) ~= 1)
                sel_source(i,j) = 1;

                %set values in mask within range to 1 --> no slides in this
                %area

                Imin=i-imask;Imax=i+imask;
                Jmin=j-jmask;Jmax=j+jmask;
                if (Imin < 1) Imin=1; end
                if (Imax > n_lat) Imax=n_lat; end
                if (Jmin < 1) Jmin=1; end
                if (Jmax > n_lon) Jmax=n_lon; end
                mask(Imin:Imax,Jmin:Jmax) = 1;
                %remove selected cell
                source(i,j) = 0;
            end 
        end %interation through colums
    end %iteration through rows

    %spread selected source cells
    temp_spreaded = climada_ls_propagation(sel_source,mult_flow,horDist,verDist,vmax,phi,iT,perWT);
    
    %%
    %reconnect spreaded slides and derive temporary area and lgt
    temp_source_ID = zeros(size(elevation));
    temp_source_ID = sel_source.*snapStart;
    [temp_IDfield,tempS] = climada_ls_slideScores(lon,lat,temp_spreaded,temp_source_ID,elevation,cell_area,field);
    
    temp_startlon = [tempS.startlon];
    temp_endlon = [tempS.endlon];
    temp_startlat = [tempS.startlat];
    temp_endlat = [tempS.endlat];
    temp_ID = [tempS.(field)];
    temp_ar = [tempS.area];
    temp_lgt = [tempS.length];
    
    temp_IDfield(temp_IDfield==0) = 100000;
    
    %find corresponding index in S of slides in sel_source; considering
    %'field' and save in array
    c = bsxfun(@eq,temp_ID,[snapS.(field)]');
    [sidx,~] = find(c); %gives field number (ID) of slide which shall be use in S.(field) to assign values
    %write lenght and area in vector
    ID(sidx) = temp_ID;
    lgt(sidx) = temp_lgt;
    area(sidx) = temp_ar;
    startlon(sidx) = temp_startlon;
    endlon(sidx) = temp_endlon;
    startlat(sidx) = temp_startlat;
    endlat(sidx) = temp_endlat;
    %merge intensity (max when overlapping) and ID (sum...)
    spreaded = max(spreaded,temp_spreaded);
    IDfield = min(IDfield,temp_IDfield);
   
    %set to zero for next round of selection
    mask(:) = 0;
    sel_source(:) = 0;
    
    %progress management
    msgstr = sprintf('\t%i of %i slides propagated...',numsource-numel(find(source==1)),numsource);
    fprintf(format_str,msgstr); % write progress to stdout
    format_str=[repmat('\b',1,length(msgstr)) '%s']; % back to begin of line
    %disp([num2str(numel(numsource)-numel(find(source==1))) filesep num2str(numel(numsource))]);
     
end %while --> ends when no sources left

fprintf(' done\n');


%write polyline (X,Y) coordinates in structure, together with ID, area and length
X = [startlon;endlon]';
Y = [startlat;endlat]';

%set coordinates of removed slides to nan
r_idx = find([snapS.removed]~=0);
X(r_idx,:) = nan;
Y(r_idx,:) = nan;
dum = [snapS.(field)];
ID(r_idx) = dum(r_idx);


for i = 1:numel(ID)
    S(i).Geometry = 'Line';
    S(i).X = X(i,:);
    S(i).Y = Y(i,:);
    S(i).removed = snapS(i).removed;
    S(i).(field) = ID(i);
    S(i).area = area(i);
    S(i).length = lgt(i);
end

%set ID to zero
IDfield(IDfield==100000) = 0;


end