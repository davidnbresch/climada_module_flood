function dist_m = climada_geo_distance_2(lonlat1,lonlat2)
% calculate distance between to points in meters
% MODULE:
%   flood
% NAME:
%   climada_geo_distance
% PURPOSE:
%   calculate distance between two points or a point and a series of points
% CALLING SEQUENCE:
%   climada_geo_distance(lonlat1, lonlat2);
% EXAMPLE:
%   climada_geo_distance(0.0,45.0,1.0,45.0); % two points
%   climada_geo_distance(0.0,45.0,[1.0 2.0 2.0],[45.0 45.0 46.0]); % point and a series of points
% INPUTS:
%   lon1,lat1: longitude and latitude of first point
%   lon2,lat2: longitude and latitude of second point
%       or, if vectors of the same length, series of points
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   dist_m: distance(s) between points in [m]
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20091227
% Lea Mueller, muellele@gmail.com, 20151106, move to flood
%-

dist_m=[]; % init

% check arguments
if ~exist('lonlat1'),fprintf('ERROR: enter lonlat1\n');return;end
if ~exist('lonlat2'),fprintf('ERROR: enter lonlat2\n');return;end
% if length(lon2)~=length(lat2),fprintf('ERROR: 2nd point vector not same length\n');return;end

% check that lonlat1 has dimension 1x2
[i1 j1] = size(lonlat1);
if i1~=1 || j1 ~=2, fprintf('ERROR: lonlat1 must have dimension 1x2 \n'), end

% check that lonlat2 has dimension Nx2
[i2 j2] = size(lonlat2);
if j2 ~=2, fprintf('ERROR: lonlat2 must have dimension Nx2 \n'), end

% fill lonlat1 into lon1 and lat1 and replicate to the dimension of lonlat2
% so that one can easily calculate with lonlat2
lon1 = repmat(lonlat1(1,1),i2,1);
lat1 = repmat(lonlat1(1,2),i2,1);

% fill lonlat2 into lon2 and lat2
lon2 = lonlat2(:,1);
lat2 = lonlat2(:,2);


% PARAMETERS
%
% constant to convert one degree latitude to km
degree2km = 111.12;

dist_m    =  sqrt( ((lon2-lon1).*cos(lat1./180.*pi) ).^2 + ...
                    (lat2-lat1)                      .^2)...
             .* 111.12 * 1000; % [in meter]
         

return
