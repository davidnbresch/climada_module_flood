function spread = climada_ls_susceptibility(lon,lat,elevation,source_areas,...
    exponent,dH,v_max,phi,friction,delta_i,perWt,d2s)

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
%   d2s:        (0/1): set to 1 if distance to source area should be
%               calculated (in [m]). Used in climada_ls_propagation where
%               the minimum distance to source area for a single slide is
%               calculated. When overlaying the all slides during one event
%               the min/max distance is taken.
%               If set 0 (default): a zeromatrix is returned.
%   single      (0/1): if set to 1 (default) --> each landslide is propagated
%               individually and the max intensity taken when all slides
%               are overlaid. Warning: calculation is slow.
%               if set 0: (not implemented at the moment) all slides are spreaded at once (faster).
%               Warning: can lead to questionable results where slides flow
%               into each other.
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
% Thomas Rölli, thomasroelli@gmail.com, 20180406, calculate distance to
%  source

%% ToDO:            Matrizen Kontrollieren ---> init at #132;  --> done
%                   speicher von resultaten (alle events --> done
%                   Distanz kontrollieren --> done
%                   Anzahl überflüsse berechnen --> done (kontrolle noch
%                   nötig)
%                   evt. implementieren dass source areas gemerged
%                   werden-->nur ein durchgang benötigt
%                   merge events --> slides müssen nur einmal berechent
%                   werden --> done (aber umsetzung nicht optimal)
%                   Gefahrenkarte erstellen --> aus allen slide und events
%                   --> partly done
%                   evt. implementieren, dass wenn flow mehrmals durchlaufen wird durchlaufen
%                   --> save --> use afterwards
%%

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
if ~exist('d2s', 'var'), d2s = []; end
%if ~exist('single', 'var'), single = []; end

%PARAMETERS
if isempty(d2s); d2s = 0; end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculation of needed data 
%%% calculate multidirectional outflow proportion
% (tan(beta_i)^x/sum(tan(beta_i)^x(i= 1 to 8))
mult_flow = climada_ls_multipleflow(lon,lat,elevation,exponent,dH);

%calculate horizontal and vertical distance to each neighbour --> needed
%when source area is spreaded downstream for energy calculation in
%climada_ls_propagation
[~,horDist,verDist] = climada_centroids_gradients(lon,lat,elevation);

%%
%%%%%%%%%%%%%%%%%%
%merging of events

sum_merged_src = sum(source_areas,3);
merged_src = sum_merged_src>0; %source areas when all events merged

source_areas = merged_src; %if each event should be calculated individually --> comment this line
n_events = numel(source_areas(1,1,:));

%%
%%%%%%%%%%%%%%%%%%%%%%%
%initiation of matrixes

%single
single_source = zeros(size(lon)); %saves source areas of single slide
single_intensity = zeros(size(lon)); %saves intensity of single slide
single_dist2source = zeros(size(lon)); %saves minimum distance to source area of single slide

%event
event_source = zeros(size(source_areas)); %stores source areas of each event
event_dist2source = zeros(size(lon)); %saves minimum distance to source of single event
event_intensity = zeros(size(lon)); %saves maximum intensity of single event
event_overflownbr = zeros(size(lon)); %saves how often a cell is affected during an event

%all
allEvents_dist2source = zeros(size(source_areas));% saves max Intensity of each event
allEvents_intensity = zeros(size(source_areas));% saves max Intensity of each event
allEvents_overflownbr = zeros(size(source_areas));% saves number of overflows of each event

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%propagation of single slides, all events and summarizing of results
cnt=0; %remove afterwards
for n = 1:n_events %iteration through events
    %set intensity and distance to zero for new event
    event_intensity = zeros(size(lon)); %stores maximum intensity of one event
    %add 100000 --> when merging minimum is taken 
    event_dist2source = zeros(size(lon))+100000; %sores minimum distance to each source area of one event
    event_source = source_areas(:,:,n); %stores source areas of each event
    for i = 1:numel(lat(:,1)) %iteration through rows
        for j = 1:numel(lon(1,:)) %iteration through collums
            if event_source(i,j)
                cnt = cnt+1;
                if (mod(cnt,100) == 0)
                    cnt 
                end
                %set single source to zero and select next source area (set
                %1) and propagate single slide 
                single_source = zeros(size(lon));
                single_source(i,j) = 1;
                [single_intensity,single_dist2source] = climada_ls_propagation(single_source,mult_flow,horDist,verDist,v_max,phi,delta_i,perWt,d2s);
                %merge single slide with event matrix (maximum of intensity, minium of
                %distance)
                event_intensity = max(event_intensity,single_intensity);
                single_dist2source(single_dist2source == 0) = 100000;
                event_dist2source = min(event_dist2source,single_dist2source);
                event_overflownbr = event_overflownbr + (single_intensity>0);
            end
        end %collums
    end %rows
    %save results from event
    allEvents_intensity(:,:,n) = event_intensity;
    %10000er values were not overflowed
    event_dist2source(event_dist2source == 100000) = 0;
    allEvents_dist2source(:,:,n) = event_dist2source;
    allEvents_overflownbr(:,:,n) = event_overflownbr;
end %events
%%



disp('test');













end
