function climada_ls_hazard_plot(lon_vec,lat_vec,field_vec,n_event)

% Plots landslide hazard event
% MODULE:
%   flood
% NAME:
%   climada_ls_hazard_plot
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
% Thomas Rölli, thomasroelli@gmail.com, 20180206, init

global climada_global
if ~climada_init_vars,return;end

% check arguments
if ~exist('lon_vec', 'var'), lon_vec = []; end
if ~exist('lat_vec', 'var'), lat_vec = []; end
if ~exist('field_vec', 'var'), field_vec = []; end
if ~exist('n_event', 'var'), n_event = 1; end

%get dimension of grid field from lon/lat coordinates
%and reshap needed vectors --> easier to handel in grid format than in
%vector; only possible for regular placed gridpoints
n_lon = numel(unique(lon_vec));
n_lat = numel(unique(lat_vec));
lon = reshape(lon_vec,n_lat,n_lon);
lat = reshape(lat_vec,n_lat,n_lon);
field = reshape(field_vec(n_event,:),n_lat,n_lon);


figure
surface(lon,lat,field)