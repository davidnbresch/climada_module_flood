function hazard = climada_hazard_encode_distance(hazard,entityORassetsORcentroids,cutoff)
% climada
% MODULE:
%   flood
% NAME:
%   climada_hazard_encode_distance
% PURPOSE:
%   Convert landslide hazard intensity (landslide depth(m))
%   into distance between an asset and the nearest centroid with nonzero
%   intensity. Distance is transformed into intensity in the form of (y = mx+b), 
%   which means increasing intensity with increasing damage, 
%   intensity = 1 - 1/cutoff * distance_m
%   Minimum distance is set to  1m, which translates into the maximum intensity 1.
%   The lat/lon coordinates of the hazard are overwritten by the asset
%   lat/lon coordinates
%   A default distance cutoff, when it can be assumed that no asset is
%   affected anymore is introduced at 1000 m.
% CALLING SEQUENCE:
%   hazard = hazard_distance_convert(hazard,entity,cutoff)
% EXAMPLE:
%   hazard = hazard_distance_convert(hazard,entity,250)
% INPUTS:
%   hazard: landslide hazard structure, intensity given as meters soildepth
%   entity: climada entity structure, with .assets.lon and .assets.lat
% OPTIONAL INPUT PARAMETERS:
%   cutoff: default 1000m, can be set to other value
% OUTPUTS:
%   hazard: a climada structure with lat/lon that equal the entity.assets.lat/lon
%           and intensity as 1-1/cutoff*distance_m (-)
% MODIFICATION HISTORY:
% Jacob Anz, j.anz@gmx.net, 20150708, initial
% Lea Mueller, muellele@gmail.com, 20150713, intensity as 1-1/cutoff*distance_m instead of distance
% Gilles Stassen, gillesstassen@hotmail.com, 20150803, faster (~30x) alternative to knnsearch; argin entity -> entityORassetsORcentroids
% Gilles Stassen, 20150805, elevation cutoff. Entity points higher in elevation than land slide excluded
% Lea Mueller, muellele@gmail.com, 20150915, bugfix when calculate distance_m from intensity
% Lea Mueller, muellele@gmail.com, 20151106, move to flood
%-

% init global variables
global climada_global
if ~climada_init_vars,return;end % init/import global variables

% poor man's version to check arguments
if ~exist('hazard',                     'var'), hazard                      = []; end
if ~exist('entityORassetsORcentroids',  'var'), entityORassetsORcentroids   = []; end
if ~exist('cutoff',                     'var'), cutoff                      = []; end

if isempty(hazard)
    hazard = climada_hazard_load; 
    if isempty(hazard),return; end
end % prompt for and load hazard, if empty

if isempty(entityORassetsORcentroids),
    entity = climada_entity_load;
    if isempty(entity),return; end
elseif isfield(entityORassetsORcentroids,'assets')      % input is entity
    entity = entityORassetsORcentroids; 
elseif isfield(entityORassetsORcentroids,'Value')       % input is assets
    entity.assets = entityORassetsORcentroids;
elseif isfield(entityORassetsORcentroids,'centroid_ID') % input is centroids
    entity.assets.lon = entityORassetsORcentroids.lon;
    entity.assets.lat = entityORassetsORcentroids.lat;
    entity.assets.elevation_m = entityORassetsORcentroids.elevation_m;
else
    cprintf([1 0 0],'ERROR: invalid input\n')
    return
end
clear entityORassetsORcentroids % prompt for and load entity, if empty

% set cutoff value, default is 1000m, all values exceeding this treshold will be set to 0
if isempty(cutoff), cutoff = 1000; end

hazard = climada_hazard2octave(hazard); % Octave compatibility for -v7.3 mat-files


%get hazard structure
hazard_distance = hazard;

% overwrite lon and lat with entity lon/lats
hazard_distance.lon        = entity.assets.lon;
hazard_distance.lat        = entity.assets.lat;
hazard_distance.intensity  = sparse(hazard.event_count, numel(entity.assets.lon));
hazard_distance.distance_m = sparse(hazard.event_count, numel(entity.assets.lon));
hazard_distance.cutoff_m   = cutoff;
hazard_distance.comment    = 'intensity, as tranformed distance, y = 1-(1/cutoff)*distance_m';
hazard_distance.units      = 'm/m';
hazard_distance.centroid_ID = 1:length(entity.assets.lon);
hazard_distance.peril_ID   = 'LS';

% init intensity matrix where we will fill in all the intensity values (y =mx+b, with x = distance_m)
intensity_matrix = zeros(hazard.event_count, numel(entity.assets.lon));

% remove fields that are not needed anymore
if isfield(hazard_distance,'source'          ), hazard_distance = rmfield(hazard_distance,'source'          ); end
if isfield(hazard_distance,'deposit'         ), hazard_distance = rmfield(hazard_distance,'deposit'         ); end
if isfield(hazard_distance,'slide_ID'        ), hazard_distance = rmfield(hazard_distance,'slide_ID'        ); end
if isfield(hazard_distance,'factor_of_safety'), hazard_distance = rmfield(hazard_distance,'factor_of_safety'); end


% sparse to full matrix
intensity_full = full(hazard.intensity);

% % find nonzero elements
% [event_indx,location_indx,s] = find(hazard.intensity);
 
stats_toolbox = 0;
% try
%     knnsearch(rand(1,2),rand(1,2)); % just to test
%     stats_toolbox = 1;
% catch
%     cprintf([ 1 0.5 0],'WARNING: no access to statistics toolbox, see line 141 in code for details on workaround\n')
%     stats_toolbox = 0;
% end

n_assets = numel(entity.assets.lon);

% init watibar
t0       = clock;
msgstr   = sprintf('encoding hazard to hazard distance');
mod_step = 10;
if climada_global.waitbar
    fprintf('%s (updating waitbar every 100th event)\n',msgstr);
    h        = waitbar(0,msgstr);
    set(h,'Name','Encoding hazard to hazard distance');
else
    fprintf('%s (waitbar suppressed)\n',msgstr);
    format_str='%s';
end
    
%loop over all events
for i=1:length(hazard_distance.event_ID)
    
    %find nonzero values
    nonzero_indx = find(intensity_full(i,:));

    if ~isempty(nonzero_indx)               
        % create array that contain longitude information in the first
        % column and latitude in the second column, works independent of
        % original lon/lat dimension
        hazard_lon_lat = [reshape(hazard.lon(nonzero_indx),numel(nonzero_indx),1) reshape(hazard.lat(nonzero_indx),numel(nonzero_indx),1)];
        entity_lon_lat = [reshape(entity.assets.lon,numel(entity.assets.lon),1) reshape(entity.assets.lat,numel(entity.assets.lat),1)];
       
        elev_check = 0;
        if isfield(entity.assets,'elevation_m') && isfield(hazard,'elevation_m')
            elev_check  = 1;
          	hazard_elev = [reshape(hazard.elevation_m(nonzero_indx),numel(nonzero_indx),1)];
            entity_elev = [reshape(entity.assets.elevation_m,numel(entity.assets.elevation_m),1)];
        end
        % find closest hazard centroid and calculate distance in meters to it
        %[indx_, distance_m] = knnsearch(hazard_lon_lat,entity_lon_lat,'Distance',@climada_geo_distance_2); 
        distance_m = [];
        if stats_toolbox 
            [~, distance_m] = knnsearch(hazard_lon_lat,entity_lon_lat,'Distance',@climada_geo_distance_2); 
        else
            % this seems to be much faster (~30x !!!)
            distance_m = zeros(n_assets,numel(nonzero_indx));
            for pt_i = 1:size(hazard_lon_lat,1)
                distance_m(:,pt_i) = climada_geo_distance_2(hazard_lon_lat(pt_i,:),entity_lon_lat);
                if elev_check,  distance_m(entity_elev>hazard_elev(pt_i),pt_i) = inf;   end
            end
            distance_m = min(distance_m,[],2);
        end
        %set minimum distance to 1m
        distance_min = 1;
        distance_m(distance_m<=distance_min) = 1;  
       
        % transform distance to an intensity in the form of (y = mx+b), which means 
        % increasing intensity with increasing damage,  
        intensity = 1 - 1./cutoff * distance_m;
        
        %apply cutoff, set all exceeding values to 0
        intensity(distance_m>=cutoff) = 0;   
        
        %new hazard with distance and lon/lat at asset location
        intensity_matrix(i,:) = intensity;
        
        % the progress management
        if mod(i,mod_step)==0
            mod_step = 100;
            t_elapsed       = etime(clock,t0)/i;
            n_remaining     = length(hazard_distance.event_ID)-i;
            t_projected_sec = t_elapsed*n_remaining;
            if climada_global.waitbar
                waitbar(i/length(hazard_distance.event_ID),h,msgstr); % update waitbar
            else
                if t_projected_sec<60
                    msgstr = sprintf('encoding hazard, est. %3.0f sec left (%i/%i events)',t_projected_sec, i,length(hazard_distance.event_ID));
                else
                    msgstr = sprintf('encoding hazard, est. %3.1f min left (%i/%i events)',t_projected_sec/60, i,length(hazard_distance.event_ID));
                end
                fprintf(format_str,msgstr); % write progress to stdout
                format_str = [repmat('\b',1,length(msgstr)) '%s']; % back to begin of line
            end
        end
    end
end
    

if climada_global.waitbar
    close(h) % dispose waitbar
else
    fprintf(format_str,''); % move carriage to begin of line
end

% write into hazard structure
hazard_distance.intensity = sparse(intensity_matrix);

% calculate distance as well (in meters)
hazard_distance.distance_m = sparse((1 - hazard_distance.intensity) .*cutoff);
% hazard_distance.distance_m = sparse(max(cutoff - hazard_distance.intensity,0));


% overwrite hazard with hazard_distance
hazard = hazard_distance;







% if plot_on
%     figure
%     scatter3(entity.assets.lon,entity.assets.lat,entity.assets.Value/max(entity.assets.Value),'.')
%     hold on
%     scatter3(hazard_distance.lon', hazard_distance.lat',hazard_distance.intensity(1,:),'*')
%     %set 0 to Nan, plots only the first event
%     testhaz=hazard.intensity(1,:);
%     testhaz(testhaz==0)=nan;
%     scatter3(hazard.lon,hazard.lat,testhaz,'.')
%     legend('assets','distance','muslide location','assets','distance');
%     sprintf('done');
% end


%%
%hazard_distance.intensity=sparse(hazard_distance.intensity);

