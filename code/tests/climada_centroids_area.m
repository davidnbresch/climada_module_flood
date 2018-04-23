function cell_area = climada_centroids_area(lon,lat,dem)

% Calculate the area of each centroid in a rectangular grid
% MODULE:
%   flood
% NAME:
%   climada_centroids_area
% PURPOSE:
%   Calculation of the area of each centroid in grid by considering the
%   slope to its 8 neighbours. Thereby, the area of each subtriangle is
%   calculated by the side-angle-side method. It is assumed that the grid
%   is rectangle and that dx (resolution in x direction [m]) and dy are the
%   same for each cell.
% CALLING SEQUENCE:
%   climada_centroids_area(lon,lat,elevation)
% EXAMPLE:
%   cell_area = climada_centroids_area(lon,lat,elevation)
% INPUTS: 
%   lat:    Latitude in gridded format (nxm) --> coordinates of DEM raster
%   lon:    longitute in gridded format (nxm) --> coordinates of DEM raster
%   dem:    digital elevation model in gridded format (nxm)
% OPTIONAL INPUT PARAMETERS:
%    
% OUTPUTS:
%    cell_area: nxm matrix of area of each centroids in rectangular grid in
%               m^2.
% MODIFICATION HISTORY:
% Thomas Rölli, thomasroelli@gmail.com, 20180423, init

global climada_global
if ~climada_init_vars, return; end

% check arguments
if ~exist('lon', 'var'), lon = []; end
if ~exist('lat', 'var'), lat = []; end
if ~exist('dem', 'var'), dem = []; end

% PARAMETERS 
if isempty(lon); return; end
if isempty(lat); return; end
if isempty(dem); return; end

%calculate vertical and horizontal distance in longitudinal/latidudinal direction and to diagonal
[~,horDist,verDist] = climada_centroids_gradients(lon,lat,dem);

dy = unique(horDist(:,:,1));
dx = unique(horDist(:,:,3));
unitarea = dx*dy;

%calculate distance of a centroid to each of its neighbours
slopDist = sqrt(horDist.^2 + verDist.^2)/2; %divided by 2 

%calculate areas of subtriangle by side-angle-side method
%area = (a*b*sin(C))/2
angle = 45; %angle between two side
area_triangle = zeros(size(slopDist));
dim = [1 2;2 3;3 4;4 5;5 6;6 7;7 8;8 1];%for selectation of corresponding sides
for c = 1:8
    area_triangle(:,:,c) = (slopDist(:,:,dim(c,1)).*slopDist(:,:,dim(c,2))*sind(angle))/2;
end

%sum up subtriangle to get area of cell
cell_area = sum(area_triangle,3);

%set boarder cells to unit area dx*dy
n_lat = numel(cell_area(:,1));
n_lon = numel(cell_area(1,:));
cell_area(1,:)=unitarea; cell_area(n_lat,:)=unitarea; cell_area(:,1)=unitarea; cell_area(:,n_lon)=unitarea;



end