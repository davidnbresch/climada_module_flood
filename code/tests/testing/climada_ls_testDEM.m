function [centroids,hazard] = climada_ls_testDEM(minmaxdegrees,resolution,...
    n_events)

% Creats a centroids and hazard set with a artificial digital elevation model and source
% areas of landslides
% test studies
% MODULE:
%   flood
% NAME:
%   climada_ls_testDEM
% PURPOSE:
%   Returns a centriods and hazard structure which a test DEM and random or
%   defined source areas of landslides
% CALLING SEQUENCE:
%   
% EXAMPLE:
%   climada_ls_multipleflow('','','',true)
% INPUTS:
%   minmaxdegrees   vector with four elemenets [minlon maxlon minlat maxlat] 
%                   defines the rectangle in which DEM shall be generated.
%   resolution:     resolution of the DEM, in arcsec
%   n_events:       number of generated events in hazard
% OPTIONAL INPUT PARAMETERS:
%
% OUTPUTS:
%   centroids:      set with lon, lat and elevation_m (artificial for
%                   tests)
%   hazard:         set with lon, lat and intensity (for n_events)
% MODIFICATION HISTORY:
% Thomas Rölli, thomasroelli@gmail.com, 20180215, init

global climada_global
if ~climada_init_vars, return; end

% check arguments
if ~exist('minmaxdegrees', 'var'), minmaxdegrees = [1 1+0.2 1 1+0.2]; end
if ~exist('resolution', 'var'), resolution = 3; end
if ~exist('n_events', 'var'), n_events = 1; end

%init centroids and hazard
centroids = [];
hazard = [];

%check if rectangle is valid
if ~(isnumeric(minmaxdegrees) && numel(minmaxdegrees) == 4)
    return
end

%generate lon lat matrixes
lon = (minmaxdegrees(1):resolution/3600:minmaxdegrees(2))';
lat = minmaxdegrees(3):resolution/3600:minmaxdegrees(4);
[lon,lat] = meshgrid(lon,lat);

%generation of elevation model
%z = 1/2.*lat;
z = 10000000000*exp(-lat.*15);

%generation of source areas
%random generiert
% source = logical(zeros(n_events,numel(lon)));
% rand = randi(10000,n_events,numel(lon));
% source(rand>=9999) = 1;

% select one source cell
source = logical(zeros(numel(lat(:,1)),numel(lon(1,:)),n_events));
source(2,round(numel(lon(1,:))/2)) = 1;

%reshape and assign data to hazard and centroids
n_ele=numel(lon);
for i = 1:n_events
    intensity(n_events,:) = reshape(source(:,:,n_events),1,n_ele);
end
%source(:,round(numel(lon(1,:))/2)) = 1;

centroids.lon = reshape(lon,1,n_ele);
centroids.lat = reshape(lat,1,n_ele);
centroids.elevation_m = reshape(z,1,n_ele);

hazard.lon = reshape(lon,1,n_ele);
hazard.lat = reshape(lat,1,n_ele);
hazard.intensity = intensity;
hazard.event_count = n_events;


end