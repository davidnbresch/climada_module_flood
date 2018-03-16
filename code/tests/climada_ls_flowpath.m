function spread = climada_ls_flowpath(lon,lat,elevation,source_areas,...
    exponent,dH,v_max,phi,friction)

% MODULE:
%   flood
% NAME:
%   climada_ls_hazard_multipleflow
% PURPOSE:
%   Process the flow path of shallow landslides with given starting
%   position (from source_areas). 
% CALLING SEQUENCE:
%   climada_ls_hazard_sets
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
%   test(removed):       If set true: a centroids structure with lon, lat and
%               elevation_m is constructed according to a test DEM (see
%               climada_ls_testDEM) (removed
% OUTPUTS:
%   mult_flow:  8-D matrix with outflow proportion in each direction (each
%               neighbour-cell)
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



%%% calculate multidirectional outflow proportion
% (tan(beta_i)^x/sum(tan(beta_i)^x(i= 1 to 8))
mult_flow = climada_ls_multipleflow(lon,lat,elevation,exponent,dH);

%calculate horizontal and vertical distance to each neighbour --> needed
%when source area is spreaded downstream 
[~,horDist,verDist] = climada_centroids_gradients(lon,lat,elevation);

%landslide path and runout distance is calculated and the intensity
%spreaded according to the multiple flow path and a simplified friction model 
spread = climada_ls_spread(source_areas,mult_flow,horDist,verDist,v_max,phi,friction);

% %make some plots just to test some things
% figure
% surf(lon,lat,elevation,spread_noFri(:,:,1));
% figure
% surf(lon,lat,elevation,spread(:,:,1));
% 
% 
% %affected area: plot to see where intensity values are greater than 0
% aff_area = zeros(n_lon,n_lat);
% aff_area(spread(:,:,1)>0) =1;
% figure
% surf(lon,lat,elevation,aff_area)
% %surf(lon,lat,elevation,aff_area,'LineStyle','none')
% 
% %affected area: plot without friction
% aff_area(spread_noFri(:,:,1)>0) =1;
% figure
% surf(lon,lat,elevation,aff_area)
% 
% disp('hier');













end
