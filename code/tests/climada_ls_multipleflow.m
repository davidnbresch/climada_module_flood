function mult_flow = climada_ls_multipleflow(lon,lat,elevation,exponent,flat_areas)

% Calculating of outflow proportion from each cell to its neighbours.
% Calculation according to Horton et al. (2013)
% MODULE:
%   flood
% NAME:
%   climada_ls_hazard_multipleflow
% PURPOSE:
%   
% CALLING SEQUENCE:
%   climada_ls_hazard_sets
% EXAMPLE:
%   
% INPUTS: 
%   centroids: a climada centroids stucture (including topographical
%              information)
% OPTIONAL INPUT PARAMETERS:
%   lon/lat:    longitudinal/latitudinal coordinates in grid format
%               -->lat(i,j)
%   elevation:  the elevation according to coordinate in grid format -->
%               elevation(i,j)
%   exponent:   variable exponent; is controlling the divergence of the flow
%               x=1: the spreading is similar to the multiple flow direction
%               x towards infinity: spreading similar to the single flow direction
%               Claessens et al. (2005) suggest x=4 for debris flow. For
%               shallow landslides 25 set by default
%   flat_areas: 0/1. if set to 1 --> searching for a path within flat
%               areas. Uses routeflats in TopoToolbox v1.06
% OUTPUTS:
%   mult_flow:  8-D matrix with outflow proportion in each direction (each
%               neighbour-cell)
%  
% MODIFICATION HISTORY:
% Thomas Rölli, thomasroelli@gmail.com, 20180201, init
% Thomas Rölli, thomasroelli@gmail.com, 20180202, calculation of outflow
%   proportion
% Thomas Rölli, thomasroelli@gmail.com, 20180305, multiple flow algorithm as
%   seperate function
% Thomas Rölli, thomasroelli@gmail.com, 20180313, flow through flat areas,
%   not working at the moment

global climada_global
if ~climada_init_vars, return; end

% check arguments
if ~exist('lon', 'var'), lon = []; end
if ~exist('lat', 'var'), lat = []; end
if ~exist('exponent', 'var'), exponent = []; end
if ~exist('flat_areas', 'var'), flat_areas = []; end


% PARAMETERS 
if isempty(exponent); exponent = 25; end
if isempty(flat_areas); flat_areas = 0; end


%calculate gradients from each cell to its 8 neighbours
[gradients,horDist,verDist] = climada_centroids_gradients(lon,lat,elevation);

% mutiplying by -1 --> outflow should be positive; neighbours with inflow are set to zero
gradients = gradients*-1;
gradients(gradients < 0) = 0; 

%%% calculate sum of all outflow cells
gradients_sum = sum(gradients.^exponent,3);
gradients_sum(gradients_sum==0) = 1; %prevent division by 0

%%% calculate multidirectional outflow proportion
% (tan(beta_i)^x/sum(tan(beta_i)^x(i= 1 to 8))
mult_flow = (gradients.^exponent)./gradients_sum;



if flat_areas
    % find path through flat areas --> see routeflats in TopoToolbox v1.06 for
    % more information
    [IXf,IXn] = routeflats(elevation,'multi');

    %get neighbour direction 1-8 from given index
    dif = IXn-IXf;
    nlat = numel(lat(:,1));
    n_ind = [+1 nlat+1 nlat nlat-1 -1 -nlat-1 -nlat -nlat+1];

    for i=1:numel(IXf)
        %get neighbour direction 1-8 from given index
        c = find(n_ind==dif(i));
        temp_mult = mult_flow(:,:,c);
        temp_mult(IXf(i))=rand(1,1);
        mult_flow(:,:,c) = temp_mult;
    end
end

end