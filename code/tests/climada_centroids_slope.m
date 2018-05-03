function slope_degree = climada_centroids_slope(lon,lat,elevation)

% calculate the neighborhood slope angle of a digital elevation model
% MODULE:
%   flood
% NAME:
%   climada_centroids_slope
% PURPOSE:
%   Returns the numerical steepest downward slope of a digital elevation
%   model using the neighbourhood slope angle (for more information see:
%   http://www.gisagmaps.com/neighborhood-slope/ and Burrough, P. A., and
%   McDonell, R. A., 1998. Principles of Geographical Information Systems 
%   (Oxford University Press, New York))
% CALLING SEQUENCE:
%   climada_centroids_slope(lon,lat,elevation)
% EXAMPLE:
%   
% INPUTS: 
%   lon/lat:    longitudinal/latitudinal coordinates in grid format
%               -->lat(n,m)
%   elevation:  the elevation according to coordinates in grid format -->
%               elevation(n,m). If the grid includes some NaNs it can lead
%               to problems
% OPTIONAL INPUT PARAMETERS:
%    
% OUTPUTS:
%    slope_degree: matrix(n,m) with slope in degrees
% MODIFICATION HISTORY:
%  Thomas Rölli, thomasroelli@gmail.com, 20180501, init

global climada_global
if ~climada_init_vars, return; end

% check arguments
if ~exist('lon', 'var'), lon = []; end
if ~exist('lat', 'var'), lat = []; end
if ~exist('elevation', 'var'), elevation = []; end

%get tan(slope) --> gradients
[~,horDist,~] = climada_centroids_gradients(lon,lat,elevation);

shift_matrix = [-1 0;-1 -1;0 -1;1 -1;1 0;1 1;0 1;-1 1];
elevation8d = zeros([size(elevation) 8]);
for c = 1:8
    elevation8d(:,:,c) = circshift(elevation,shift_matrix(c,:));
end

%calculate east-west gradient dz/dx
x_operator = reshape([0 1 2 1 0 -1 -2 -1],[1 1 8]);
y_operator = reshape([2 1 0 -1 -2 -1 0 1],[1 1 8]);

dx = bsxfun(@times,elevation8d,x_operator);
dx_dz = sum(dx,3)./(8*horDist(:,:,3)); %matrix 3 --> horizontal distance towards the East

dy = bsxfun(@times,elevation8d,y_operator);
dy_dz = sum(dy,3)./(8*horDist(:,:,1)); %matrix 1 --> horizontal distance towards the North

slope_gradient = sqrt(dx_dz.^2 + dy_dz.^2);

slope_degree = atand(slope_gradient);

% set slopes of boarder cells to zero
slope_degree([1 numel(lat(:,1))],:) = 0; %upper boarder
slope_degree(:,[1 numel(lon(1,:))]) = 0; %side boarder

end
