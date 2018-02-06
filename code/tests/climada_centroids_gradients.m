function gradients = climada_centroids_gradients(lon,lat,elevation)

% Calculation of gradients of each cell to its corresponding 8 neighbours.
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
%   elevation: DEM data with the elevation in meters in grid-format
%   lon: longitudinal information (position of centroids and therefore DEM grids) 
%   lat: latitudinal information (position of centroids and therefore DEM grids)
% OPTIONAL INPUT PARAMETERS:
%   
% OUTPUTS:
%   gradients:  8-D matrix with gradients in each direction (towards each
%               neighbour-cell)
%  
% MODIFICATION HISTORY:
% Thomas Rölli, thomasroelli@gmail.com, 20180206, init

global climada_global
if ~climada_init_vars, return; end

% check arguments
if ~exist('elevation', 'var'), return; end

% PARAMETERS
%for calculations
deg_km = 111.12; %length of 1 degree on Earth


%%%% calculate tan(beta_i)%%%%
%slope of all center-cell to all 8 neigbour-cells (direction i)

%calculate distance in longitudinal/latidudinal direction and to diagonal
%neighbour-cells
%it is asumed that the longitudinal distance doesn't change with latitude
%(mean latitude is taken)
dy = min(diff(lat(:,1)))*(deg_km * 1000);
dx = min(diff(lon(1,:)))*cosd(mean(lat(:,1)))*(deg_km * 1000); 
dxdy = sqrt(dx^2+dy^2);

%shif matrix such that neigbour cell is shifted on the centre cell (with circshift)
%starting at 12 o'clock and proceeding clockwise
shift_matrix = [1 0;1 -1;0 -1;-1 -1;-1 0;-1 1;0 1;1 1];
hor_dist = [dy,dxdy,dx,dxdy,dy,dxdy,dx,dxdy];

% gradients(:,:,1) corresponds to Northerly (12o'clock) cell, then going clockwise
gradients = zeros(numel(lat(:,1)),numel(lat(1,:)),8);

for c = 1:8
    gradients(:,:,c) = (circshift(elevation,shift_matrix(c,:))-elevation)...
        /hor_dist(c);
end

% set gradients which are directed out of the study region to zero; at
% boarder cells
gradients(1,:,[1 2 8]) = 0; %upper boarder
gradients(numel(lat(:,1)),:,[4 5 6]) = 0; %lower boarder
gradients(:,1,[6 7 8]) = 0; %left boarder
gradients(:,numel(lon(1,:)),[2 3 4]) = 0; %right boarder

