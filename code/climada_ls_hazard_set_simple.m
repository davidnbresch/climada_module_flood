function hazard  = climada_ls_hazard_set_simple(centroids,n_events,hazard_set_file,wiggle_factor,TWI_condition,wiggle_factor_slope,slope_condition,...
                   n_downstream_cells,focus_area,polygon_correction,random_trigger_condition)
% Generate a very simple landslide hazard set
% MODULE:
%   flood
% NAME:
%   climada_ls_hazard_set_simple
% PURPOSE:
%   Generate land slide hazard set based on topographical information 
%   (slope and topographical wetness index)
% PREVIOUS STEP:
%   centroids_fl_prepare
% CALLING SEQUENCE:
%   hazard = climada_ls_hazard_set_simple(centroids, n_events, wiggle_factor,TWI_condition,wiggle_factor_slope,slope_condition)
% EXAMPLE:
%   hazard = climada_ls_hazard_set_simple(centroids)
% INPUTS:
%   centroids:  a climada centroids stucture including topographical
%   information
% OPTIONAL INPUT PARAMETERS:
%   n_events: number of events
%   hazard_set_file: the name (and path, optional) of the hazard set file
%       If no path provided, default path ../data/hazards is used (and name
%       can be without extension .mat)
%       > promted for if not given
%   wiggle_factor: an array, default is 0.35, to modify topographical
%   wetness factor, which is a number between 0 and 1.4
%   TWI_condition: an array, default is 0.95, to define where a landslide occurs
%   wiggle_factor_slope: an array, default is 0.35, to modify the slope
%         factor, which is a number between 0 and 1
%   slope_condition: an array, default is 0.45, to define where a landslide occurs
%   n_downstream_cells: number of downstream cells where the landslide occurs as well
%   focus_area: a polygon to define the focus area (with focus_area.lon, focus_area.lat),
%         landslides only in the given area will be filtered
%   polygon_correction: a polygon to define an area where less landlislides should occur
%   random_trigger_condition: a number between 0 and 1, 1 prevents all
%         landslide in the polygon_correction area, 0 does not inhibit any
%         landlides in the polygon_correction area
% OUTPUTS:
%   hazard: a climada hazard structure
%       peril_ID: 'LS'
%       date: the creation date of the set
%       intensity(event_i,centroid_i),sparse: the hazard intensity of event_i at
%           centroid_i
%       frequency(event_i): the frequency of each event
%       matrix_density: the density of the sparse array hazard.intensity
%       filename: the filename of the hazard event set (if passed as a
%           struct, this is often useful)
% MODIFICATION HISTORY:
% Lea Mueller, muellele@gmail.com, 20150911, init
% Lea Mueller, muellele@gmail.com, 20150920, add polygon_correction with a random trigger condition, to reduce landslides in a given polygon
% -

hazard = []; % init

global climada_global
if ~climada_init_vars, return; end

% check arguments
if ~exist('centroids', 'var'), centroids = []; end
if ~exist('n_events', 'var'), n_events = []; end
if ~exist('hazard_set_file', 'var'), hazard_set_file = []; end
if ~exist('wiggle_factor', 'var'), wiggle_factor = []; end
if ~exist('TWI_condition', 'var'), TWI_condition = []; end
if ~exist('wiggle_factor_slope', 'var'), wiggle_factor_slope = []; end
if ~exist('slope_condition', 'var'), slope_condition = []; end
if ~exist('n_downstream_cells', 'var'), n_downstream_cells = []; end
if ~exist('focus_area', 'var'), focus_area = []; end
if ~exist('polygon_correction', 'var'), polygon_correction = []; end
if ~exist('random_trigger_condition', 'var'), random_trigger_condition = []; end


if isempty(centroids) 
    return
end

if isempty(n_events); n_events = 100; end
if isempty(wiggle_factor); wiggle_factor = 0.35; end
if isempty(TWI_condition); TWI_condition = 0.95; end
if isempty(wiggle_factor_slope); wiggle_factor_slope = 0.2; end
if isempty(slope_condition); slope_condition = 0.45; end
if isempty(n_downstream_cells); n_downstream_cells = 5; end
if isempty(random_trigger_condition); random_trigger_condition = 0; end


% prompt for hazard_set_file if not given
if isempty(hazard_set_file) % local GUI
    hazard_set_file      = [climada_global.data_dir filesep 'hazards' filesep 'LSXX_hazard.mat'];
    [filename, pathname] = uiputfile(hazard_set_file, 'Save LS hazard set as:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        hazard_set_file = fullfile(pathname,filename);
    end
end

% complete path, if missing
[fP,fN,fE]=fileparts(hazard_set_file);
if isempty(fP),hazard_set_file=[climada_global.data_dir filesep 'hazards' filesep fN fE];end


%% load template for hazard-structure
hazard_example_file = [climada_global.data_dir filesep 'hazards' filesep 'TCNA_today_small.mat'];
if exist(hazard_example_file,'file')
    load(hazard_example_file)
else
    fprintf('No hazard example found to be loaded. \n')
end
hazard_ex = hazard;

% overwrite template hazrd with actual information
hazard.lon = centroids.lon;
hazard.lat = centroids.lat;
hazard.centroid_ID = 1:numel(hazard.lon);
hazard.peril_ID         = 'LS';
hazard.orig_years       = 10000;
hazard.orig_event_count = n_events;
hazard.event_count      = n_events;
hazard.event_ID         = 1:n_events;
hazard.orig_event_flag  = ones(1,n_events);
hazard.yyyy             = ones(1,n_events);
hazard.mm               = ones(1,n_events);
hazard.dd               = ones(1,n_events);
hazard.intensity        = sparse(n_events,numel(hazard.lon));
hazard.name             = cell(1,n_events);
hazard.frequency        = ones(1,n_events)/hazard.orig_years;
hazard.comment          = centroids.comment;
hazard.date             = datestr(now);
hazard.units            = 'binary';
hazard.orig_yearset     = [];
hazard.filename         = hazard_set_file;
hazard.matrix_density   =  0.01; % estimate

n_centroids = length(centroids.centroid_ID);

if ~isfield(centroids,'slope_deg')
    fprintf('Add topographical characteristics to the centroids, based on elevation\n')
    centroids = centroids_TWI(centroids);
end

% calculate slope_factor as cos(slope)/sin(slope)
slope_factor = 1./(cosd(centroids.slope_deg) ./ sind(centroids.slope_deg));
slope_factor(isinf(slope_factor)) = 0;
slope_factor(slope_factor>0.55) = 0.6;
if ~isfield(centroids,'slope_factor')
    centroids.slope_factor = slope_factor;
end

% normalize TWI
if ~isfield(centroids,'TWI_norm')
    % TWI_norm = centroids.TWI/10;
    % TWI_norm(TWI_norm>0.85) = 0.85;
    TWI_norm = centroids.TWI/10;
    TWI_norm(isnan(TWI_norm)) = 0;
    centroids.TWI_norm = TWI_norm;
end

% % slope factor
% titlestr = 'Slope factor (-)';
% miv = 0;
% mav = 0.6;
% fig = climada_figuresize(0.5,0.6);
% plotclr(centroids.lon, centroids.lat, centroids.slope_factor, marker,markersize,cbar_on,miv,mav);
% title(titlestr); axis(axlim); box on; climada_figure_scale_add('',7,1)
% pdf_filename = sprintf('LS_slope_factor.pdf');
% print(fig,'-dpdf',[ls_dir pdf_filename])
% 
% % TWI norm
% titlestr = 'TWI normalised';
% miv = 0;
% mav = 1;
% fig = climada_figuresize(0.5,0.6);
% plotclr(centroids.lon, centroids.lat, centroids.TWI_norm, marker,markersize,cbar_on,miv,mav);
% title(titlestr); axis(axlim); box on; climada_figure_scale_add('',7,1)
% pdf_filename = sprintf('LS_TWI_normalised.pdf');
% print(fig,'-dpdf',[ls_dir pdf_filename])

if ~isempty(polygon_correction)
    lon_lat_polygon = climada_concatenate_lon_lat(polygon_correction.lon,polygon_correction.lat);
    lon_lat = climada_concatenate_lon_lat(centroids.lon,centroids.lat);
    needs_correction = inpoly(lon_lat,lon_lat_polygon);
    %sum(needs_correction)
    needs_correction = find(needs_correction);
end

% create TWI wiggle as the sum of TWI_norm + TWI_delta
% create slope wiggle as the sum of slope_factor + slope_delta
% TWI_wiggle = zeros(size(hazard.intensity)); %init
% slope_wiggle = zeros(size(hazard.intensity)); %init

TWI_delta = rand(size(hazard.intensity)) * wiggle_factor;
slope_delta = rand(size(hazard.intensity)) * wiggle_factor_slope;

for e_i = 1:n_events
    TWI_wiggle = centroids.TWI_norm + TWI_delta(e_i,:);
    slope_wiggle = centroids.slope_factor + slope_delta(e_i,:);
        
    % check where landslides occur
    ls_occurence(e_i,:) = TWI_wiggle>TWI_condition & slope_wiggle>slope_condition ;
    
    if ~isempty(polygon_correction)
        random_trigger = rand(numel(needs_correction),1);
        random_trigger(random_trigger>random_trigger_condition) = 1;
        random_trigger(random_trigger<=random_trigger_condition) = 0;
        not_triggered = needs_correction(~random_trigger);
        ls_occurence(e_i,not_triggered) = 0;
    end
    
end


% expand landslides to following n_cells downstream
if n_downstream_cells>1
    if ~isfield(centroids,'sink_ID_10')
        centroids = climada_flow_find(centroids);
    end
    if n_downstream_cells>10
        n_downstream_cells = 10;
        fprintf('Maximum numbers of downstream cells is 10.\n')
    end
    
    msgstr   = sprintf('Expand landslide to %d downstream cells for %i events ... ',n_downstream_cells,n_events);
    mod_step = 10; % first time estimate after 10 assets, then every 100
    if climada_global.waitbar
        fprintf('%s (updating waitbar with estimation of time remaining every 50th event)\n',msgstr);
        h        = waitbar(0,msgstr);
        set(h,'Name','Expand landslide to ownstream cells');
    else
        fprintf('%s (waitbar suppressed)\n',msgstr);
        format_str='%s';
    end
    
    % loop over all events
    for e_i = 1:n_events
        is_event = ls_occurence(e_i,:);
        centroid_list = find(is_event);
        
        % loop over centroids, that are sliding
        for i = centroid_list
            selected_sinks = centroids.sink_ID_10(i,:);
            % take only a given number of downstream cells
            if numel(selected_sinks)>n_downstream_cells
                 selected_sinks = selected_sinks(1:n_downstream_cells);
            end
            if ~isempty(selected_sinks)
                is_sink = []; %init
                for s_i = 1:numel(selected_sinks)
                    if selected_sinks(s_i)>0
                        is_sink(s_i) = find(selected_sinks(s_i) == centroids.centroid_ID);
                    end
                end
                %is_sink = ismember(centroids.centroid_ID,selected_sinks);
                ls_occurence(e_i,is_sink) = ones(1,numel(is_sink));
            end
        end
        
        % the progress management
        if mod(e_i,mod_step)==0
            mod_step          = 50;
            msgstr = sprintf('%i/%i events',e_i,n_events);
            if climada_global.waitbar
                waitbar(e_i/n_events,h,msgstr); % update waitbar
            else
                fprintf(format_str,msgstr); % write progress to stdout
                format_str=[repmat('\b',1,length(msgstr)) '%s']; % back to begin of line
            end
        end
    
    end
    ls_occurence = logical(ls_occurence);
end

hazard.intensity = sparse(ls_occurence);

if ~isempty(focus_area)
    fprintf('Filter out landslides in the given focus area.\n')
    if isfield(focus_area,'lon') & isfield(focus_area,'lat')
        polygon = climada_concatenate_lon_lat(focus_area(1).lon, focus_area(1).lat);
        
    elseif isfield(focus_area,'X') & isfield(focus_area,'Y')
        polygon = climada_concatenate_lon_lat(focus_area(1).X, focus_area(1).Y);
        
    elseif isnumeric(focus_area) % it is already formatted as a polygon
        [i, j] = size(focus_area)
        if j == 2 & i>2
            polygon = focus_area;
        end
    else
        fprintf('Please check the input of the focus area.\n')
        return
    end
    hazard_lon_lat = climada_concatenate_lon_lat(hazard.lon, hazard.lat);
    is_inside = inpoly(hazard_lon_lat,polygon);
    if any(is_inside)
        hazard.lon = hazard.lon(is_inside);
        hazard.lat = hazard.lat(is_inside);
        hazard.centroid_ID = 1:numel(hazard.lon);
        hazard.intensity = hazard.intensity(:,is_inside);
    end
end


fprintf('Save landslide (LS) hazard set as %s\n',hazard_set_file);
save(hazard_set_file,'hazard')




