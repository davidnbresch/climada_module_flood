function [hazard,centroids] = climada_ls_hazard_sets(centroids,srtm1,n_events,set_files,...
    wiggle_factor_TWI,condition_TWI,wiggle_factor_slope,condition_slope,spread_exponent,...
    v_max,phi_friction)

% Generate a landslide hazard set.
% MODULE:
%   flood
% NAME:
%   climada_ls_hazard_set
% PURPOSE:
%   Test version for generating shallow landslide hazards (some parts are
%   adapted from climada_ls_hazard_set.m
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
%   [hazard, centroids] = function [hazard,centroids] = climada_ls_hazard_test(centroids,n_events,...
%                                 wiggle_factor_TWI,condition_TWI,wiggle_factor_slope,condition_slope)

% EXAMPLE:
%   %%TEST, for a region around Sarnen in Switzerland (Kt. Obwalden):
%   [hazard,centroids]=climada_ls_hazard_sets([8.2456-.05 8.2456+.05 46.8961-.05 46.8961+.05],0,100,'_LS_Sarnen');
%   with srtm1-data
%   [hazard,centroids]=climada_ls_hazard_sets([8.2456-.05 8.2456+.05 46.8961-.05 46.8961+.05],1,100,'_LS_Sarnen_srtm1');
% INPUTS:
%   centroids:  a climada centroids stucture (ideally including topographical
%       information) or a rectangle to define lon/lat box, if not given, the
%       user can select a rectangle by first selecting a country and then drawing
%       a rectangle on a map
% OPTIONAL INPUT PARAMETERS:
%   n_events: number of events
%   set_file: the name (and path, optional) of the hazard set and centroids
%       files. If no path provided, default path ../data/hazards is used (and name
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
%   srtm1:                if set to 1 srtm1-data is used to get the
%                         elevation of the terrain. 0 by default: srtm3 data is used 
%   spread_exponent       variable exponent; is controlling the divergence of the flow
%                         x=1: the spreading is similar to the multiple flow direction
%                         x towards infinity: spreading similar to the single flow direction
%                         Claessens et al. (2005) suggest x=4 for debris flow. For
%                         shallow landslides 25 set by default
%   v_max:                describes maximal possible velocity of slide. If velocity is 
%                         exceeded v_max is taken. Should keep energy amounts within reasonal values
%                         and therefore prevent improbalbe runout distances
%   phi_friction:         angle of reach, angle of the line connnecting the source area to
%                         the most distant point reached by the slide, along its path.
%                         Factor controlls maximum possible runout distance
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
% Thomas Rölli, thomasroelli@gmail.com, 20180201, set up test programm 

% init
hazard = []; fig = [];

global climada_global
if ~climada_init_vars, return; end

% check arguments
if ~exist('centroids', 'var'), centroids = []; end
if ~exist('srtm1', 'var'), srtm1 = 0; end
if ~exist('n_events', 'var'), n_events = []; end
if ~exist('set_files', 'var'), set_files = []; end
if ~exist('wiggle_factor', 'var'), wiggle_factor_TWI = []; end
if ~exist('TWI_condition', 'var'), condition_TWI = []; end
if ~exist('wiggle_factor_slope', 'var'), wiggle_factor_slope = []; end
if ~exist('slope_condition', 'var'), condition_slope = []; end
if ~exist('n_downstream_cells', 'var'), n_downstream_cells = []; end
if ~exist('focus_area', 'var'), focus_area = []; end
if ~exist('polygon_correction', 'var'), polygon_correction = []; end
if ~exist('random_trigger_condition', 'var'), random_trigger_condition = []; end
if ~exist('check_plot', 'var'), check_plot = 0; end
if ~exist('spread_exponent', 'var'), spread_exponent = []; end
if ~exist('v_max', 'var'), v_max = []; end
if ~exist('phi_friction', 'var'), phi_friction = []; end

% prompt for set_files if not given
if isempty(set_files) % local GUI
    hazard_set_file      = [climada_global.data_dir filesep 'hazards' filesep 'LSXX_hazard.mat'];
    centroids_set_file    = [climada_global.data_dir filesep 'centroids' filesep 'LSXX_centroids.mat'];
    %prompt for hazard file
    [filename,pathname] = uiputfile(hazard_set_file, 'Save LS hazard set as:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        hazard_set_file = fullfile(pathname,filename);
    end
    %prompt for centroids file
    [filename,pathname] = uiputfile(centroids_set_file, 'Save LS centroids set as:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        centroids_set_file = fullfile(pathname,filename);
    end
else
    hazard_set_file      = [climada_global.data_dir filesep 'hazards' filesep set_files '_hazard.mat'];
    centroids_set_file   = [climada_global.data_dir filesep 'centroids' filesep set_files '_centroids.mat'];
end

%create centroids and add elevation information to it (from SRTM 90)
if isnumeric(centroids) && numel(centroids) == 4
    % we have a box that defines where centroids should be created on 90m
    % resolution (given by SRTM)
    % return
    centroids = climada_centroids_elevation_add('',centroids,0,srtm1);
end

if isempty(centroids) 
    % create centroids by asking user for a country and to define a
    % rectangle region on the figure
    % return
    centroids = climada_centroids_elevation_add('','',0,srtm1);
end

%calculate TWI, aspect, slope
centroids = climada_centroids_TWI_calc_v2(centroids);

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


%create hazard set file and assess susceptibility of shallow landslides --> get trigger areas for e.g.
%100 events (1/0)
hazard = climada_ls_hazard_trigger(centroids,n_events,...
    wiggle_factor_TWI,condition_TWI,wiggle_factor_slope,condition_slope);

%get dimension of grid field from lon/lat coordinates
%and reshap needed vectors --> easier to handel in grid format than in
%vector; and needed in climada_ls_flowpath
n_lon = numel(unique(centroids.lon));
n_lat = numel(unique(centroids.lat));
lon = reshape(centroids.lon,n_lat,n_lon);
lat = reshape(centroids.lat,n_lat,n_lon);
elevation = reshape(centroids.elevation_m,n_lat,n_lon);
intensity = logical(zeros(n_lat,n_lon,hazard.event_count));
for i = 1:hazard.event_count
    intensity(:,:,i) = reshape(hazard.intensity(i,:),n_lat,n_lon);
end

%assess flow path of landslide; spread source areas in intensity donwstream
spread = climada_ls_flowpath(lon,lat,elevation,intensity,spread_exponent,v_max,phi_friction);



if isempty(hazard),return;end % Cancel pressed in climada_ls_hazard_set_binary or failed

% save hazard and centroids
hazard.filename = hazard_set_file;
centroids.filename = centroids_set_file;
fprintf('Save landslide (LS) hazard and centroids set as %s\n and %s\n',hazard_set_file,centroids_set_file);
save(hazard_set_file,'hazard')
save(centroids_set_file,'centroids')

end