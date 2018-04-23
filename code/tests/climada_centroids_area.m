function climada_centroids_area(lon,lat,elevation)

% Script to calibrate flow path parameters
% MODULE:
%   flood
% NAME:
%   climada_ls_FPcalibration
% PURPOSE:
%   
% CALLING SEQUENCE:
%   
% EXAMPLE:
%   
% INPUTS: 
%     
% OPTIONAL INPUT PARAMETERS:
%    
% OUTPUTS:
%    
% MODIFICATION HISTORY:
% Thomas Rölli, thomasroelli@gmail.com, 20180423, init

%remove afterwards
load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_srtm1_centroids.mat','centroids')
n_lon = numel(unique(centroids.lon));
n_lat = numel(unique(centroids.lat));
lon = reshape(centroids.lon,n_lat,n_lon);
lat = reshape(centroids.lat,n_lat,n_lon);
elevation = reshape(centroids.elevation_m,n_lat,n_lon);

global climada_global
if ~climada_init_vars, return; end

% check arguments
if ~exist('lon', 'var'), lon = []; end
if ~exist('lat', 'var'), lat = []; end
if ~exist('elevation', 'var'), elevation = []; end

% PARAMETERS 
if isempty(lon); return; end
if isempty(lat); return; end
if isempty(elevation); return; end

%calculate vertical and horizontal distance in longitudinal/latidudinal direction and to diagonal
[~,horDist,verDist] = climada_centroids_gradients(lon,lat,elevation);

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
cell_area(1,:)=unitarea; cell_area(n_lat,:)=unitarea; cell_area(:,1)=unitarea; cell_area(:,n_lon)=unitarea;



end