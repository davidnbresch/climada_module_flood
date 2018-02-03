function mult_flow = climada_ls_multipleflow(centroids,exponent)

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
%   exponent:   variable exponent; is controlling the divergence of the flow
%               x=1: the spreading is similar to the multiple flow direction
%               x towards infinity: spreading similar to the single flow direction
%               Claessens et al. (2005) suggest x=4 for debris flow 
%  
% OUTPUTS:
%   mult_flow:  8-D matrix with outflow proportion in each direction (each
%               neighbour-cell)
%  
% MODIFICATION HISTORY:
% Thomas Rölli, thomasroelli@gmail.com, 20180201, init
% Thomas Rölli, thomasroelli@gmail.com, 20180202, calculation of outflow
%   proportion

%remove afterwards; load centroids
load('C:\Users\Simon Rölli\Desktop\climada\climada_data\centroids\_LS_Sarnen_centroids.mat')

global climada_global
if ~climada_init_vars, return; end

% check arguments
if ~exist('centroids', 'var'), centroids = []; end
if ~exist('exponent', 'var'), exponent = []; end

% PARAMETERS 
if isempty(exponent); exponent = 4; end
%for calculations
deg_km = 111.12; %length of 1 degree on Earth

%extract needed data from centroids and reshap in grid fromat --> easier to
%handle; only possible for regular placed gridpoints
n_lon = numel(unique(centroids.lon));
n_lat = numel(unique(centroids.lat));
lon = reshape(centroids.lon,n_lat,n_lon);
lat = reshape(centroids.lat,n_lat,n_lon);
elevation = reshape(centroids.elevation_m,n_lat,n_lon);


%%%% calculate tan(beta_i)%%%%
%slope of all center-cell to all 8 neigbour-cells (direction i)

%calculate distance in longitudinal/latidudinal direction and to diagonal
%neighbour-cells
%it is asumed that the longitudinal distance doesn't change with latitude
%(mean latitude is taken)
dy = min(diff(unique(centroids.lat)))*(deg_km * 1000);
dx = min(diff(unique(centroids.lon)))*cosd(mean(centroids.lat))*(deg_km * 1000); 
dxdy = sqrt(dx^2+dy^2);

%shif matrix such that neigbour cell is shifted on the centre cell (with circshift)
%starting at 12 o'clock and proceeding clockwise
shift_matrix = [1 0;1 -1;0 -1;-1 -1;-1 0;-1 1;0 1;1 1];
hor_dist = [dy,dxdy,dx,dxdy,dy,dxdy,dx,dxdy];

% gradients(:,:,1) corresponds to Northerly (12o'clock) cell, then going clockwise
gradients = zeros(n_lat,n_lon,8);

for c = 1:8
    gradients(:,:,c) = (circshift(elevation,shift_matrix(c,:))-elevation)...
        /hor_dist(c);
end

% mutiplying by 1 --> outflow should be positive; neighbours with inflow are set to zero
gradients = gradients*-1;
gradients(gradients < 0) = 0; 


%%% calculate sum of all outflow cells
gradients_sum = sum(gradients,3);
gradients_sum(gradients_sum==0) = 1; %prevent division by 0


%%% calculate multidirectional outflow proportion
% (tan(beta_i)^x/sum(tan(beta_i)(i= 1 to 8))^x

mult_flow = (gradients)./(gradients_sum);
%mult_flow = (gradients.^exponent)/(gradients_sum.^exponent);











end
