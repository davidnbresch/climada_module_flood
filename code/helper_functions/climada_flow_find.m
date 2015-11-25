function centroids = climada_flow_find(centroids)
% Find flow direction for every centroid 
% MODULE:
%   flood
% NAME:
%   climada_flow_find
% PURPOSE:
%   Find flow direction for every centroid, find next 10 centroids, based
%   on sink_ID
% PREVIOUS STEP:
%   centroids = climada_centroids_TWI_calc(centroids)
% CALLING SEQUENCE:
%   centroids = climada_flow_find(centroids)
% EXAMPLE:
%   centroids = climada_flow_find(centroids)
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
%   centroids: a climada centroids structure, including TWI properties,
%   especially .sink_ID
% OUTPUTS:
%   centorids: centroid with field .sink_ID_10 defining the 10 next centroids
% MODIFICATION HISTORY:
% Lea Mueller, muellele@gmail.com, 20150920, init
% Lea Mueller, muellele@gmail.com, 20151125, rename to climada_centroids_TWI_calc from centroids_TWI
% -


global climada_global
if ~climada_init_vars, return; end

n_centroids = numel(centroids.sink_ID);
n_sinks = 10;
sink_ID = zeros(n_centroids,n_sinks); % init

msgstr   = sprintf('Find flow direction for %i centroids ... ',n_centroids);
mod_step = 10; % first time estimate after 10 assets, then every 100
if climada_global.waitbar
    fprintf('%s (updating waitbar with estimation of time remaining every 100th centroid)\n',msgstr);
    h        = waitbar(0,msgstr);
    set(h,'Name','Find flow direction for centroids');
else
    fprintf('%s (waitbar suppressed)\n',msgstr);
    format_str='%s';
end


for c_i = 1:n_centroids
    next_sink_ID = centroids.sink_ID(c_i);
    
    for c_ii = 1:n_sinks
        if ~isempty(next_sink_ID) & ~isnan(next_sink_ID)         
            sink_ID(c_i,c_ii) = next_sink_ID;
            next_sink_ID = centroids.sink_ID(sink_ID(c_i,c_ii) == centroids.centroid_ID);
            
            %this_elevation = centroids.elevation_m(sink_ID(c_i,c_ii) == centroids.centroid_ID);
            %next_elevation = centroids.elevation_m(centroids.centroid_ID == next_sink_ID);
            
            % %local minimum, do not use
            % if this_elevation<next_elevation 
            %     % go back to last sink ID
            %     previous_sink_ID = sink_ID(c_i,c_ii-1);
            % 
            %     % find all IDs that flow into this local minimum
            %     have_similar_sink = centroids.sink_ID == sink_ID(c_i,c_ii);
            %     centroids.centroid_ID(have_similar_sink)
            % 
            %     have_similar_sink_2 = ismember(centroids.sink_ID,centroids.centroid_ID(have_similar_sink));
            %     centroids.centroid_ID(have_similar_sink_2)
            % end
        end
    end
    
    % the progress management
    if mod(c_i,mod_step)==0
        mod_step          = 100;
        msgstr = sprintf('%i/%i centroids',c_i,n_centroids);
        if climada_global.waitbar
            waitbar(c_i/n_centroids,h,msgstr); % update waitbar
        else
            fprintf(format_str,msgstr); % write progress to stdout
            format_str=[repmat('\b',1,length(msgstr)) '%s']; % back to begin of line
        end
    end
    
end
    
centroids.sink_ID_10 = sink_ID;


if climada_global.waitbar
    close(h) % dispose waitbar
else
    fprintf(format_str,''); % move carriage to begin of line
end


% figure
% is_sink = ismember(centroids.centroid_ID,sink_ID(1,:));
% plot(centroids.lon(is_sink),centroids.lat(is_sink),'o-')
