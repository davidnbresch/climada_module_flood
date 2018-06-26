function [modS,obsS] = climada_ls_probAss(lon,lat,elevation,subS,snapS,...
    edges_vmax,edges_phi,count_phi,num_slides,fig)

% Generates probabilistic parameter sets (vmax, phi) from
% climada_ls_probVmaxPhi, propagates n randomly chosen slides in snapS
% and compares it with the landslide inventory in subS (area and length)
% MODULE:
%  
% NAME:
%  climada_probSet_ls
% PURPOSE:
%  Can be used to assess the performace of the flow path algorithm for a 
%  probabilistic parameter set (different vmax phi values). Probabilistic
%  set shall represent the natural variablity and uncertainty in the flow
%  path parameters. A good parameter sets leads to a good representation of
%  the overall distribution, range of the area and length of the shallow
%  landslides from the inventory
% CALLING SEQUENCE:
%   
% EXAMPLE:
%   
% INPUTS:
%  lon,lat,elevation: DEM grid data. Same DEM as used to produce snapS.
%  subS:    Shape structure with 'Line' as Geometriy. Derived from original polygon shapes
%           with climada_calibrat_orgS2raster
%           Need to include fields: X,Y,length,area,reachAngle, andremoved (1/0) --> slides 
%           which shall not be considered
%  snapS:   Shape structure with 'Line' as Geometry. Derived from subS
%           after snapping to nearest grid cell in coarser grid (with
%           climada_calibration_snap). Need to include fields: X and Y
% OPTIONAL INPUT PARAMETERS:
%   edges_vmax: vector with 4 elements. First and last element gives range
%               of possible vmax values. From the first to the second edge
%               in the vector, the probability is linear increasing.
%               Between element 2 and 3, the probability is constant and
%               between 3 and 4, linear decreasing. Outside of the range,
%               (element 1 and 4) the probabiliy is zero. Default = [1 4 8 11]
%   edges_phi:  Vector which was derived from [~,edges_phi] =
%               histcounts(...). (use Defines the edges of the histogram of
%               angle of reach (phi) of landslide inventory. if not
%               provided, histcounts is applied on [subS.reachAngle] with a
%               barwidth of 2.5.
%   count_phi:  Vector which was derived from [count_phi,~] =
%               histcounts(...). Values should be normalized (use
%               'Normalization', 'probability'). Defines bar height of
%               histogram of angle of reach (phi) of landslide inventory.
%   num_slides: number of slides which shall be randomly choosen from snapS
%               and propagated with random parameter (vmax,phi).
%   fig:        plotting (1/0) of parameter distribution in
%               climada_ls_probVmaxPhi
% OUTPUTS:   
%   
% MODIFICATION HISTORY:
% Thomas Rölli, thomasroelli@gmail.com, 20180611, init


global climada_global
if ~climada_init_vars, return; end

%check arguments
if ~exist('lon') lon = []; end
if ~exist('lat') lat = []; end
if ~exist('elevation') elevation = []; end
if ~exist('subS') subS = []; end
if ~exist('snapS') snapS = []; end
if ~exist('edges_vmax') edges_vmax = []; end
if ~exist('edges_phi') edges_phi = []; end
if ~exist('count_phi') count_phi = []; end
if ~exist('num_slides') num_slides = []; end
if ~exist('fig') fig = []; end

if isempty(lon); return; end
if isempty(lat); return; end
if isempty(elevation); return; end
if isempty(subS); return; end
if isempty(snapS); return; end
if isempty(edges_vmax); edges_vmax = [1 4 8 11]; end
if numel(edges_vmax) ~= 4; return; end

reachAngle = [subS.reachAngle];
rmv = [subS.removed];
if isempty(edges_vmax) || isempty(count_phi) 
    [count_phi,edges_phi] = histcounts(reachAngle,'BinWidth',2.5,'Normalization', 'probability');
end

if isempty(num_slides); num_slides = 1000; end
if isempty(fig); fig = 0; end



%generate probabilistic set of parameters (vmax,phi)
[vmax,phi] = climada_ls_probVmaxPhi(count_phi,edges_phi,edges_vmax,10^7,'lin',fig);

%define other flow path parameters (which are kept constant)
exp = 25;
PT = 0.0003; 
perWt = [1 0.8 0.4 0 0 0 0.4 0.8];


%get data for slide propagation
elevation = fillsinks(elevation);
cell_area = climada_centroids_area(lon,lat,elevation);
%mult_flow = climada_ls_multipleflow(lon,lat,elevation);
mult_flow = climada_ls_multipleflow(lon,lat,elevation,exp);
[~,horDist,verDist] = climada_centroids_gradients(lon,lat,elevation);


%create vector which chose randomly num_slides times a slide in subS and
%also randomly assigns a parameter set
slide_idx = randi([1 numel(subS)],1,num_slides);
para_idx = randi([1 numel(vmax)],1,num_slides);

vmax = vmax(para_idx);
phi = phi(para_idx);

format_str = '%s';

%propagate all slides with corresponding parameter set
for i=1:num_slides
    sl_idx = slide_idx(i);
    if rmv(sl_idx)~=1 %just for slides which are not removed
        start = zeros(size(lat));
        %extract coordinates from polyline structure
        x = snapS(sl_idx).X;
        y = snapS(sl_idx).Y;
        %find indices for start cell
        [~,idxJ] = ismember(x(1),lon(1,:));
        [~,idxI] = ismember(y(1),lat(:,1));
        start(idxI,idxJ) = 1;
        
        %propagate slide
        spreaded = climada_ls_propagation(start,mult_flow,horDist,verDist,vmax(i),phi(i),PT,perWt);
        %calculate length area 
        source_ID = start*snapS(sl_idx).ID;
        [~,tempS] = climada_ls_slideScores(lon,lat,spreaded,source_ID,elevation,cell_area);
        
        %write data in modS and obsS
        X = [tempS.startlon tempS.endlon];
        Y = [tempS.startlat tempS.endlat];
        modS(i).ID = tempS.ID; modS(i).X = X; modS(i).Y = Y; modS(i).length = tempS.length; modS(i).area = tempS.area;
        modS(i).vmax = vmax(i); modS(i).phi = phi(i);
        
        obsS(i) = subS(sl_idx);
        modS(i).removed = 0;
        
        
        if mod(i,100) == 0
            msgstr = sprintf('\t%i of %i slides propagated...',i,num_slides);
            fprintf(format_str,msgstr); % write progress to stdout
            format_str=[repmat('\b',1,length(msgstr)) '%s']; % back to begin of line
        end
    else
        modS(i).removed = 1;
        obsS(i) = subS(sl_idx);
    end
end

%remove slides in rmv from structure --> structure can be saved in shape file
%modS(rmv).IDmodS(rmv_idx) = [];
rmv_idx = find([modS.removed]~=0);
modS(rmv_idx) = [];
obsS(rmv_idx) = [];

%add geometry field
[modS(:).Geometry] = deal('Line');

end