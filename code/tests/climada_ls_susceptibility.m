function spread = climada_ls_susceptibility(lon,lat,elevation,source_areas,...
    exponent,dH,v_max,phi,friction,delta_i,perWt)

% MODULE:
%   flood
% NAME:
%   climada_ls_susceptibility
% PURPOSE:
%   Process the flow path of shallow landslides with given starting
%   position (starting from source_areas). Each source area is propagated
%   individually. After all slides are processed, the maximum of intensity is taken when
%   several slides flew over on cell. This is repeated for each event.
%   
%   While propagating (in climada_ls_propagation), the minimum flow-distance to the
%   source area is saved and the maximum taken when including all slides
%   of an event. 
%   
%   Last but not least, number of times each cell was affected by a slide is computed (for each event). 
% CALLING SEQUENCE:
%   climada_ls_susceptibility(lon,lat,elevation,source_areas,exponent,dH,v_max,phi,friction,delta_i,perWt)
% EXAMPLE:
%   
% INPUTS: 
%   lon/lat:    longitudinal/latitudinal coordinates in grid format
%               -->lat(i,j)
%   elevation:  the elevation for coordinate points in grid format -->
%               elevation(i,j)
%   source_area: starting points of shallow landslides for serveral events
%                in grid fromat -->source_areas(i,j,event_n) = 0/1
% OPTIONAL INPUT PARAMETERS:
%   exponent:   variable exponent; is controlling the divergence of the flow
%               x=1: the spreading is similar to the multiple flow direction
%               x towards infinity: spreading similar to the single flow direction
%               Claessens et al. (2005) suggest x=4 for debris flow
%   dH:         changing the height of the central cell by a factor dH.
%               Allows smoothing of DEM and leads to a more consistent
%               spreading. It can also be used to transfer the flow
%               throught flat areas. Used in climada_centroids_gradients
%   v_max:      describes maximal possible velocity of slide. If velocity is 
%               exceeded v_max is taken. Should keep energy amounts within reasonal values
%               and therefore prevent improbalbe runout distances
%   phi:        angle of reach, angle of the line connnecting the source area to
%               the most distant point reached by the slide, along its path.
%               Factor controlls maximum possible runout distance
%   friction:   1/0 include friction/no friction while spreading
%   delta_i:    Minimum threshold which prevent propagation when deceeded. Small values
%               are removed and spreaded to the other remaining, lower
%               situated neighbouring cells 
%   perWt:      row vector with 8 elements. gives weight to each direction.
%               The persistence function aims at
%               reproducing the behaviour of inertia --> weights the flow
%               direction according to change of direction to neighbour.
%               First element represent weight to neighbour in same direction
%               as flow (0 degree), second element weights right neighbour 45
%               degrees (angle between previous direction and direction from
%               the central cell to corresponding neighbour) .. last element
%               represents weight of neighbour to the left of flow direction
%               (45 degree anticlockwise). 
% OUTPUTS:
%   sucept:     matrix (lon lat) with inte
%  
% MODIFICATION HISTORY:
% Thomas Rölli, thomasroelli@gmail.com, 20180201, init
% Thomas Rölli, thomasroelli@gmail.com, 20180202, calculation of outflow
%   proportion
% Thomas Rölli, thomasroelli@gmail.com, 20180208, implementation of flow
%   path
% Thomas Rölli, thomasroelli@gmail.com, 20180214, implementaiton of
%   friction, outflow distance
% Thomas Rölli, thomasroelli@gmail.com, 20180227, do not flip lat anymore
% Thomas Rölli, thomasroelli@gmail.com, 20180306, remove test-DEM and
%  lat/lon, elevation and intensity in grid is now demanded in gridded
%  format.
% Thomas Rölli, thomasroelli@gmail.com, 20180403, add delta_i and perWt and
%  uses spread_v2
% Thomas Rölli, thomasroelli@gmail.com, 20180404, renamed from
%  climada_ls_flowpath to climada_ls_susceptibility

%remove afterwards; load centroids and hazard
%load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards\_LS_Sarnen_hazard.mat')
%load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards\_LS_Sarnen_srtm1_hazard.mat')

global climada_global
if ~climada_init_vars, return; end
% check arguments
if ~exist('lon', 'var'), lon = []; end
if ~exist('lat', 'var'), lat = []; end
if ~exist('elevation', 'var'), elevation = []; end
if ~exist('source_areas', 'var'), source_areas = []; end
if ~exist('exponent', 'var'), exponent = []; end
if ~exist('dH', 'var'), dH = []; end
if ~exist('v_max', 'var'), v_max = []; end
if ~exist('phi', 'var'), phi = []; end
if ~exist('friction', 'var'), friction = []; end
if ~exist('delta_i', 'var'), delta_i = []; end
if ~exist('perWt', 'var'), friction = []; end

n_events = numel(source_areas(1,1,:));

%%% calculate multidirectional outflow proportion
% (tan(beta_i)^x/sum(tan(beta_i)^x(i= 1 to 8))
mult_flow = climada_ls_multipleflow(lon,lat,elevation,exponent,dH);

%calculate horizontal and vertical distance to each neighbour --> needed
%when source area is spreaded downstream 
[~,horDist,verDist] = climada_centroids_gradients(lon,lat,elevation);

%landslide path and runout distance is calculated and the intensity
%spreaded according to the multiple flow path and a simplified friction model 
spread = climada_ls_spread(source_areas,mult_flow,horDist,verDist,v_max,phi,friction);

%with spread_v2 --> each slide is spreaded seperately

%spreaded_v2 = climada_ls_propagation(source_areas,mult_flow,horDist,verDist,v_max,phi,delta_i,perWt);
count=1;
for n = 1:n_events %iteration through events
    event_total = zeros(size(source_areas));
    for i = 1:numel(lat(:,1)) %iteration through rows
        for j = 1:numel(lon(1,:)) %iteration through collums
            if source_areas(i,j)
                %count
                count = count+1;
                single_spread = zeros(size(source_areas));
                single_spread(i,j) = 1;
                single_spreaded = climada_ls_propagation(single_spread,mult_flow,horDist,verDist,v_max,phi,delta_i,perWt);
                event_total = max(event_total,single_spreaded);
            end
        end %collums
    end %rows
    
end %events

disp('test');













end
