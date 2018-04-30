function [gradients,horDist,verDist] = climada_centroids_gradients(lon,lat,elevation,dH,lat_dep)

% Calculation of gradients of each cell to its corresponding 8 neighbours. 
% Make sure to use regular gridded data.
% MODULE:
%   flood
% NAME:
%   climada_ls_hazard_multipleflow
% PURPOSE:
%   returns tan(slope); to receive the slope in degrees use
%   atand(gradients)
% CALLING SEQUENCE:
%   climada_ls_hazard_sets
% EXAMPLE:
%   
% INPUTS: 
%   elevation: DEM data with the elevation in meters in grid-format
%   lon: longitudinal information (position of centroids and therefore DEM grids) 
%   lat: latitudinal information (position of centroids and therefore DEM grids)
% OPTIONAL INPUT PARAMETERS:
%   dH:  elevation central cell is raised by dH (in units of dem data).
%        Allows smoothing, especially for DEMs with high resolution.
%        Default: dH=0
%   lat_dep:  0/1. Default: lat_dep=0; If set to 1: The change of the
%             longitudional distance of 1 degree with latitude is
%             considered when calculating the horizontal distance between
%             the grid points.
% OUTPUTS:
%   gradients:  8-D matrix with gradients (tan(slope)) in each direction (towards each
%               neighbour-cell), gradients(:,:,1) gradient toward northern
%               neighbour, then continuing clockwise
%   horDist:    8-D Matrix with horizontal distance in each direction
%   verDist:    8-D Matrix with vertical distance in each direction
%   dH:         changing the height of the central cell by dH.
%               Allows smoothing of DEM and leads to a more consistent
%               spreading. It can also be used to transfer the flow
%               throught flat areas. Used in climada_centroids_gradients
%  
% MODIFICATION HISTORY:
% Thomas Rölli, thomasroelli@gmail.com, 20180206, init
% Thomas Rölli, thomasroelli@gmail.com, 20180227, flipped data not needed
%  anymore
% Thomas Rölli, thomasroelli@gmail.com, 20180315, add dH
% Thomas Rölli, thomasroelli@gmail.com, 20180424, add latitude-dependence
% Thomas Rölli, thomasroelli@gmail.com, 20180426, set boarder to zero

global climada_global
if ~climada_init_vars, return; end

% check arguments
if ~exist('elevation', 'var'), return; end
if ~exist('lat', 'var'), return; end
if ~exist('lon', 'var'), return; end
if ~exist('dH', 'var'), dH = []; end
if ~exist('lat_dep', 'var'), lat_dep = []; end

% PARAMETERS
%for calculations
if isempty(dH), dH = 0; end
if isempty(lat_dep), lat_dep = 0; end

%E = wgs84Ellipsoid;
% E.SemimajorAxis;
%deg_km_lat = E.SemiminorAxis*2*pi/(360*1000);
%deg_km_lon = E.SemimajorAxis*2*pi/(360*1000);

deg_km_lon = 111.32; %length of 1 degree on Earth in longitudinal direction
deg_km_lat = 111.32;
dlat = abs(min(diff(lat(:,1))));
dlon = abs(min(diff(lon(1,:))));

% gradients(:,:,1) corresponds to Northerly (12o'clock) cell, then going clockwise
gradients = zeros([size(lat) 8]);
horDist = zeros([size(lat) 8]);
verDist = zeros([size(lat) 8]);

%%
%calculate vertical distance to neighbour

%shif matrix such that neigbour cell is shifted on the centre cell (with circshift)
%starting at 12 o'clock and proceeding clockwise (12 o'clock --> toward
%North)
shift_matrix = [-1 0;-1 -1;0 -1;1 -1;1 0;1 1;0 1;-1 1];
for c = 1:8
    verDist(:,:,c) = circshift(elevation,shift_matrix(c,:))-(elevation+dH);
end
%%
%calculate horizontal distance to neighbour

%without considering lat dependency
if ~lat_dep
    %calculate distance in longitudinal/latidudinal direction and to diagonal
    %neighbour-cells
    %it is asumed that the longitudinal distance doesn't change with latitude
    %(mean latitude is taken)
    dy = dlat*(deg_km_lat * 1000)*ones(size(lat));
    dx = dlon*cosd(mean(lat(:,1)))*(deg_km_lon * 1000)*ones(size(lat)); 
    dxdy = sqrt(dx.^2+dy.^2);
    %horizontal distance
    horDist = cat(3,dy,dxdy,dx,dxdy,dy,dxdy,dx,dxdy);
elseif lat_dep
    %calculate distance in longitudinal/latidudinal direction and to diagonal
    %neighbour-cells (considering longitudinal dependency)
    dx = dlon*cosd(lat)*(deg_km_lon * 1000);
    dy = dlat*(deg_km_lat * 1000)*ones(size(lat));
    dxdy = sqrt(dx.^2+dy.^2);
    %horizontal distance
    horDist = cat(3,dy,dxdy,dx,dxdy,dy,dxdy,dx,dxdy);
end
%%
% calculate tan(beta_i)%%
%slope of all center-cell to all 8 neigbour-cells (direction i)
gradients = verDist./horDist;

% set gradients which are directed out of the study region to zero; at
% boarder cells
gradients(numel(lat(:,1)),:,[1 2 8]) = 0; %upper boarder
gradients(1,:,[4 5 6]) = 0; %lower boarder
gradients(:,1,[6 7 8]) = 0; %left boarder
gradients(:,numel(lon(1,:)),[2 3 4]) = 0; %right boarder
%same for verDist
verDist(numel(lat(:,1)),:,[1 2 8]) = 0; %upper boarder
verDist(1,:,[4 5 6]) = 0; %lower boarder
verDist(:,1,[6 7 8]) = 0; %left boarder
verDist(:,numel(lon(1,:)),[2 3 4]) = 0; %right boarder

