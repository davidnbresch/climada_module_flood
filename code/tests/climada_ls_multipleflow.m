function mult_flow = climada_ls_multipleflow(lon,lat,elevation,exponent,dH,flat_areas)

% Calculating of outflow proportion from each cell to its neighbours.
% Calculation according to Horton et al. (2013)
% 
% MODULE:
%   flood
% NAME:
%   climada_ls_hazard_multipleflow
% PURPOSE:
%   Calculates the outflow proportion of each cell to its neighbours. It
%   uses the Holmgren's multiple flow algorithm (see Horton et al. (2013): 
%   Flow-R, a model for susceptibility mapping of gravitational hazards.) B
%   using a TopoToolbox function the algorithm handels also flat areas but
%   it is not recommended to use this option (flat_areas = 1) because it is
%   not working properly. At the moment it seems to make no big difference
%   if you ignore flat areas or not, because shallow landslide normally do
%   not slide over flat areas. Nevertheless for other application (e.g.
%   calculation of flowaccumulation) consider to use the TopoToolbox
%   functions (see FLOWobj(DEM,'multi') --> uses multiple flow algorithm 
%   according to Freeman 1991) which treats flat areas properly. Before
%   using the multipleflow algorithm the DEM should be free of sinks. Use
%   the function fillsinks(DEM) in TopoToolbox to get rid of sinks.
% CALLING SEQUENCE:
%   mult_flow =
%   climada_ls_multipleflow(lon,lat,elevation,exponent,flat_areas
% EXAMPLE:
%   climada_ls_multipleflow(lon,lat,elevation,25,1)
% INPUTS: 
%   lon/lat:    longitudinal/latitudinal coordinates in grid format
%               -->lat(i,j)
%   elevation:  the elevation according to coordinates in grid format -->
%               elevation(i,j). If the grid includes some NaNs it can lead
%               to problems
% OPTIONAL INPUT PARAMETERS:
%   exponent:   variable exponent; is controlling the divergence of the flow
%               x=1: the spreading is similar to the multiple flow direction
%               x towards infinity: spreading similar to the single flow direction
%               Claessens et al. (2005) suggest x=4 for debris flow. For
%               shallow landslides: 25 set by default.
%   dH:         changing the height of the central cell by a factor dH.
%               Allows smoothing of DEM and leads to a more consistent
%               spreading. It can also be used to transfer the flow
%               throught flat areas but leads to too high values in and 
%               after flat areas. Used in climada_centroids_gradients
%   flat_areas: 0/1. if set to 1 --> searching for a path within flat
%               areas. Uses routeflats in TopoToolbox v1.06, it is
%               recommended to use this option.
% OUTPUTS:
%   mult_flow:  8-D matrix with outflow proportion in each direction (each
%               neighbour-cell). mult_flow(:,:,1) shows the outflow to
%               northern neighbour, from there the matrix (:,:,2)...
%               proceeds clockwise
%  
% MODIFICATION HISTORY:
% Thomas Rölli, thomasroelli@gmail.com, 20180201, init
% Thomas Rölli, thomasroelli@gmail.com, 20180202, calculation of outflow
%   proportion
% Thomas Rölli, thomasroelli@gmail.com, 20180305, multiple flow algorithm as
%   seperate function
% Thomas Rölli, thomasroelli@gmail.com, 20180313, flow through flat areas,
%   not working at the moment
% Thomas Rölli, thomasroelli@gmail.com, 20150315, add dH, fix problem with
%   TopoToolbox function
global climada_global
if ~climada_init_vars, return; end

% check arguments
if ~exist('lon', 'var'), lon = []; end
if ~exist('lat', 'var'), lat = []; end
if ~exist('elevation', 'var'), elevation = []; end
if ~exist('exponent', 'var'), exponent = []; end
if ~exist('dH', 'var'), dH = []; end
if ~exist('flat_areas', 'var'), flat_areas = []; end


% PARAMETERS 
if isempty(exponent); exponent = 25; end
if isempty(flat_areas); flat_areas = 0; end


%calculate gradients from each cell to its 8 neighbours
[gradients,horDist,verDist] = climada_centroids_gradients(lon,lat,elevation,dH);

% mutiplying by -1 --> outflow should be positive; neighbours with inflow are set to zero
gradients = gradients*-1;
gradients(gradients < 0) = 0; 

%%% calculate sum of all outflow cells
gradients_sum = sum(gradients.^exponent,3);
gradients_sum(gradients_sum==0) = 1; %prevent division by 0

%%% calculate multidirectional outflow proportion
% (tan(beta_i)^x/sum(tan(beta_i)^x(i= 1 to 8))
mult_flow = (gradients.^exponent)./gradients_sum;



% find path through flat areas --> see routeflats in TopoToolbox v1.06 for
% more information
if flat_areas
    %IXf includeds the indices of all the identified cells within a flat
    %area, while IXn indcludes the next neighbour where the corresponding
    %IXf should flow in order to reach the outlet of the flat.
    [IXf,IXn] = routeflats(elevation,'multi');

    %difference of two indices describes which neighbour respectively
    % matrix (e.g. (:,:,1)) the slide should flow to reach outlet
    dif = IXn-IXf;
    nlat = numel(lat(:,1));
    n_ind = [+1 nlat+1 nlat nlat-1 -1 -nlat-1 -nlat -nlat+1];

    for i=1:numel(IXf)
        %get neighbour direction 1-8 from given index
        c = find(n_ind==dif(i));
        %outflow value in corresponding cell (flat area cell)
        temp_mult = mult_flow(:,:,c);
        temp_mult(IXf(i))=1;
   
        mult_flow(:,:,c) = temp_mult;
    end
    %in flat areas: sum over neighbours can be greater than 1 (because in routflats()
    % multiple flow option was used) --> normalise again
    mult_flow_sum = sum(mult_flow,3);
    mult_flow_sum(mult_flow_sum==0) = 1;
    mult_flow = mult_flow./mult_flow_sum;
end



end