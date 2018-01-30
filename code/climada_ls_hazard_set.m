function [hazard,centroids,fig] = climada_ls_hazard_set(centroids,n_events,hazard_set_file,...
    wiggle_factor_TWI,condition_TWI, wiggle_factor_slope,condition_slope,...
    n_downstream_cells,focus_area,polygon_correction,random_trigger_condition,...
    check_plot)
% Generate a landslide hazard set.
% MODULE:
%   flood
% NAME:
%   climada_ls_hazard_set
% PURPOSE:
%   This is a all-in-one function to generate a landslide hazard set. The
%   landslide events are based solely on topographical information
%   (SRTM 90m digital elevation data) and calculates slope and
%   topographical wetness index.

%   We create first a binary hazard, encode to distance and plot figures if required.
%   invokes
%       - climada_ls_hazard_set_binary
%       - climada_hazard_encode_distance
%   and for plots
%       - climada_ls_hazard_binary_plot
%       - climada_map_plot
%       - climada_hazard_stats
%       - climada_hazard_stats_figure
%
%   Useful TEST data for the TEST area in EXAMPLE below:
%   https://www.gis-daten.ch/map/ow_natgef_planungszoneHWschutz_gwr?extentMinX=645100&extentMinY=180800&extentMaxX=679500&extentMaxY=201800&view=Ereigniskataster
%
%   previous call: e.g.create centroids, select a rectangle box to create the landslide hazard
%   centroids
%   centroids = climada_centroids_elevation_add;
% CALLING SEQUENCE:
%   [hazard, centroids, fig] = climada_ls_hazard_set(centroids,n_events,hazard_set_file,...
%                                         wiggle_factor_TWI,condition_TWI, wiggle_factor_slope,condition_slope,...
%                                         n_downstream_cells,focus_area,polygon_correction,random_trigger_condition,...
%                                         check_plot)
% EXAMPLE:
%   %%TEST, for a region around Sarnen in Switzerland (Kt. Obwalden):
%   [hazard,centroids]=climada_ls_hazard_set([8.2456-.05 8.2456+.05 46.8961-.05 46.8961+.05],100,'_LS_Sarnen_binary');
%   %[hazard,centroids]=climada_ls_hazard_set([8.2456-.05 8.2456+.05 46.8961-.05 46.8961+.05],100,'_LS_Sarnen_binary','','','','','','','','',1); % with check plots
%   checksum=sum(sum(hazard.intensity)) should be = 1.0060e+06 % check-sum (plus/minus, since radom generator)
%   figure;climada_hazard_plot_nogrid(hazard,0,2); % show max hazard intensity at each centroid
%   climada_hazard_stats(hazard); % show hazard intensity return periods
%   figure;plot3(centroids.lon,centroids.lat,centroids.elevation_m,'.r') % show terrain
%   % for a tile in San Salvador (where the code has been first developed for:
%   [hazard, centroids]  = climada_ls_hazard_set([-89.145 -89.1 13.692 13.727],'','','','','','','','','','',1)
% INPUTS:
%   centroids:  a climada centroids stucture (ideally including topographical
%       information) or a rectangle to define lon/lat box, if not given, the
%       user can select a rectangle by first selecting a country and then drawing
%       a rectangle on a map
% OPTIONAL INPUT PARAMETERS:
%   n_events: number of events
%   hazard_set_file: the name (and path, optional) of the hazard set file
%       If no path provided, default path ../data/hazards is used (and name
%       can be without extension .mat)
%       > promted for if not given
%   wiggle_factor_TWI:    an array, default is 0.35, to modify topographical
%                         wetness factor, which is a number between 0 and 1.4
%   condition_TWI:        an array, default is 0.95, to define a minimum
%                         topographical wetness index, where a landslide occurs
%   wiggle_factor_slope:  an array, default is 0.35, to modify the slope factor,
%                         which is a number between 0 and 1
%   condition_slope:      an array, default is 0.45, to define a minimum slope
%                         where a landslide occurs
%   n_downstream_cells:   number of downstream cells where the landslide is extended
%   focus_area:           a polygon to define the focus area (with focus_area.lon, focus_area.lat),
%                         landslides only in the given area will be filtered
%                         and the other areas are cut out
%   polygon_correction:   a polygon to define an area where less landlislides should occur
%   random_trigger_condition: a number between 0 and 1, 1 prevents all landslide
%                         in the polygon_correction area, 0 does not inhibit any
%                         landlides in the polygon_correction area
%   check_plot:          set to 1 if you want to see maps (elevation, slope, hazard, etc)
% OUTPUTS:
%   hazard: a climada hazard structure with binary landslide information
%       .peril_ID: 'LS'
%       .date: the creation date of the set
%       .intensity(event_i,centroid_i),sparse: the hazard intensity of event_i at
%           centroid_i, binary, 1 indicates a landslide, 0 no landslide
%       .frequency(event_i): the frequency of each event
%       .matrix_density: the density of the sparse array hazard.intensity
%       .filename: the filename of the hazard event set (if passed as a
%           struct, this is often useful)
%   centroids: a climada centroids structure with fields
%       .lon
%       .lat
%       .elevation_m, elevation in meters
%       .slope_deg, slope in degree
%       .TWI, topographical wetness index
%       .aspect_deg, spect in degree
%   fig: handle of figures (maps with elevation, slope, landslide
%        hazard,etc)
% MODIFICATION HISTORY:
% Lea Mueller, muellele@gmail.com, 20151124, init
% David N. Bresch, david.bresch@gmail.com, 20171116, made more stable, example added
% David N. Bresch, david.bresch@gmail.com, 20171117, check plots fixed
% -

% init
hazard = []; fig = [];

global climada_global
if ~climada_init_vars, return; end

% check arguments
if ~exist('centroids', 'var'), centroids = []; end
if ~exist('n_events', 'var'), n_events = []; end
if ~exist('hazard_set_file', 'var'), hazard_set_file = []; end
if ~exist('wiggle_factor', 'var'), wiggle_factor_TWI = []; end
if ~exist('TWI_condition', 'var'), condition_TWI = []; end
if ~exist('wiggle_factor_slope', 'var'), wiggle_factor_slope = []; end
if ~exist('slope_condition', 'var'), condition_slope = []; end
if ~exist('n_downstream_cells', 'var'), n_downstream_cells = []; end
if ~exist('focus_area', 'var'), focus_area = []; end
if ~exist('polygon_correction', 'var'), polygon_correction = []; end
if ~exist('random_trigger_condition', 'var'), random_trigger_condition = []; end
if ~exist('check_plot', 'var'), check_plot = 0; end


% create binary landslide hazard
[hazard_binary, centroids]  = climada_ls_hazard_set_binary(centroids,n_events,hazard_set_file,...
    wiggle_factor_TWI,condition_TWI, wiggle_factor_slope,condition_slope,...
    n_downstream_cells,focus_area,polygon_correction,random_trigger_condition);

if isempty(hazard_binary),return;end % Cancel pressed in climada_ls_hazard_set_binary or failed

% encode binary hazard to distance, so that we have the distance to
% landslides as intensity
cutoff = 1000;
hazard = climada_hazard_encode_distance(hazard_binary,centroids,cutoff);

% save hazard
[pathname, filename, f_ext] = fileparts(hazard.filename);
if ~exist(pathname,'dir')
    hazard_set_file = [climada_global.data_dir filesep 'hazards' filesep 'LSXX_hazard.mat'];
    [filename, pathname] = uiputfile(hazard_set_file, 'Save LS hazard set as:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        hazard_set_file = fullfile(pathname,filename);
    end
else
    filename = [strrep(filename,'binary','') 'distance'];
    hazard_set_file = fullfile(pathname,filename,f_ext);
end
fprintf('Save landslide (LS) hazard set (encoded to distance) as %s\n',hazard_set_file);
save(hazard_set_file,'hazard')

if check_plot
    % create landslide binary event map
    fig = climada_ls_hazard_binary_plot(hazard_binary);
    
    figure; % new figure
    % plot centroids with characteristics (elevation, slope, twi, etc)
    fieldname_to_plot = {'elevation_m' 'slope_deg' 'TWI' 'aspect_deg'};
    plot_method = 'plotclr';
    [~, fig_temp] = climada_map_plot(centroids,fieldname_to_plot,plot_method);
    
    % plot hazard statistics
    hazard.orig_years = 1000; % we set the number of years to 1000 to have nice images, but please check if this is suitable
    hazard.frequency = ones(size(hazard.event_ID))*(1./hazard.orig_years);
    
    figure; % new figure
    return_periods = [10 25 50 100 150 200];
    hazard_stats = climada_hazard_stats(hazard,return_periods,0);
    fig_temp_2 = climada_hazard_stats_figure(hazard_stats,return_periods);
    
    % concatenate figure handles
    fig = [fig fig_temp fig_temp_2];
end

end % climada_ls_hazard_set
