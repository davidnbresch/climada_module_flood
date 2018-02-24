function [gradients,horDist,verDist] = climada_centroids_gradients(lon,lat,elevation)

% Calculation of gradients of each cell to its corresponding 8 neighbours.
% Pay attention: ordered such that latitude starts from high values and
% decreases with each row (consider to use flipud()). Also make sure you
% use regular gridded data.
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
%       should start from high lat and decreases with each row
% OPTIONAL INPUT PARAMETERS:
%   
% OUTPUTS:
%   gradients:  8-D matrix with gradients in each direction (towards each
%               neighbour-cell), gradients(:,:,1) gradient toward northern
%               neighbour, then continuing clockwise
%   horDist:    8-D Matrix with horizontal distance in each direction
%   verDist:    8-D Matrix with vertical distance in each direction
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
dy = abs(min(diff(lat(:,1)))*(deg_km * 1000));
dx = abs(min(diff(lon(1,:)))*cosd(mean(lat(:,1)))*(deg_km * 1000)); 
dxdy = sqrt(dx^2+dy^2);

%shif matrix such that neigbour cell is shifted on the centre cell (with circshift)
%starting at 12 o'clock and proceeding clockwise
shift_matrix = [1 0;1 -1;0 -1;-1 -1;-1 0;-1 1;0 1;1 1];
hor_dist = [dy,dxdy,dx,dxdy,dy,dxdy,dx,dxdy];

% gradients(:,:,1) corresponds to Northerly (12o'clock) cell, then going clockwise
gradients = zeros(numel(lat(:,1)),numel(lat(1,:)),8);
horDist = zeros(numel(lat(:,1)),numel(lat(1,:)),8);
verDist = zeros(numel(lat(:,1)),numel(lat(1,:)),8);

for c = 1:8
    horDist(:,:,c) = hor_dist(c);
    verDist(:,:,c) = circshift(elevation,shift_matrix(c,:))-elevation;
    %gradients(:,:,c) = (circshift(elevation,shift_matrix(c,:))-elevation)...
    %    /hor_dist(c);
end
gradients = verDist./horDist;

% set gradients which are directed out of the study region to zero; at
% boarder cells
gradients(1,:,[1 2 8]) = 0; %upper boarder
gradients(numel(lat(:,1)),:,[4 5 6]) = 0; %lower boarder
gradients(:,1,[6 7 8]) = 0; %left boarder
gradients(:,numel(lon(1,:)),[2 3 4]) = 0; %right boarder

