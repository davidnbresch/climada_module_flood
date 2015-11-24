function fig = climada_ls_hazard_binary_plot(hazard)
% climada plot ls binary hazard
% MODULE:
%   flood
% NAME:
%   climada_ls_hazard_binary_plot
% PURPOSE:
%   This is a function that is usually called from climada_ls_hazard_set.
%   Plot binary landslide hazard, all events on one map
% PREVIOUS STEP:
%   hazard_binary = climada_ls_hazard_set_binary(centroids)
% NEXT STEP:
%   hazard = climada_hazard_encode_distance(hazard,centroids,cutoff);
% CALLING SEQUENCE:
%   fig = climada_ls_hazard_binary_plot(hazard)
% EXAMPLE:
%   fig = climada_ls_hazard_binary_plot;
% INPUTS:
%   hazard:  a climada hazard stucture with binary landslide
%   information (.intensity is either 1 or 0)
% OPTIONAL INPUT PARAMETERS:
%   none
% OUTPUTS:
%   fig: handle of map with landslide events
% MODIFICATION HISTORY:
% Lea Mueller, muellele@gmail.com, 20151124, init
% -

global climada_global
if ~climada_init_vars, return; end

% check arguments
if ~exist('hazard', 'var'), hazard = []; end

if isempty(hazard), hazard = climada_hazard_load; end   
         

% visualize binary ls hazard on a map
% ------------------------------
% delete nans if there are any
hazard.intensity(isnan(hazard.intensity)) = 0;
axlim = [min(hazard.lon) max(hazard.lon) min(hazard.lat) max(hazard.lat)];
n_events = hazard.event_count;
n_colors = jet(n_events);

% show a maximum number of events
event_step = 1;
max_events_to_show = 100;
if n_events>max_events_to_show, event_step = round(n_events/max_events_to_show);end

fig = climada_figuresize(0.5,0.7);
% plot(entity.assets.lon, entity.assets.lat,'.','linewidth',0.2,'markersize',0.8,'color',[255 64 64 ]/255);
hold on
legendstr = []; h = []; counter = 0;
for e_i = 1:event_step:n_events
    is_event = logical(hazard.intensity(e_i,:));
    if any(is_event)
        counter = counter+1;
        %hold on; plot3(hazard.lon(is_event), hazard.lat(is_event), ones(sum(is_event))*3000, 'dr','linewidth',2,'markersize',5,'color',[255 64 64 ]/255)
        h(counter) = plot(hazard.lon(is_event), hazard.lat(is_event),'dr','linewidth',2,'markersize',5,'color',n_colors(e_i,:));
        hold on; 
        %plot(polygon_canas.X, polygon_canas.Y, 'b-');
        legendstr{counter} = sprintf('Event %d',e_i);
    end
end
titlestr = 'LS hazard binary';
try titlestr = sprintf('%s hazard (%s), %d events, %d',hazard.peril_ID, hazard.units, n_events, hazard.reference_year);end
title(titlestr); axis(axlim); box on; 
climada_figure_axis_limits_equal_for_lat_lon(axlim); climada_figure_scale_add('',1,1)
legend(h,legendstr,'location','eastoutside')







