function centroids = centroids_basin_ID(centroids, res, basin_shapefile, check_plots)
% assign basin ID to centroids
% MODULE:
%   tbd
% NAME:
%	country_basin_ID
% PURPOSE:
%   Assign basin IDs to centroids based on HydroSHEDS basin outline
%   shapefiles For more information and a technical documentation see
%   http://hydrosheds.org/page/hydrobasins
% CALLING SEQUENCE:
%   centroids_basin_ID = centroids_basin_ID(basin_shapefile,centroids, check_plots)
% EXAMPLE:
%   centroids_basin_ID = centroids_basin_ID('',centroids,0)
% INPUTS:
%   centroids: Climada centroids struct; the following fields are required:
%         .lat:           Latitude
%         .lon:           Longitude
%         .centroid_ID:   centroid ID
%         .admin0_name    Country name
% OPTIONAL INPUT PARAMETERS:
%   basin_shapefile: name (including full path) of the basin shapefile to
%   be used. If left empty, the code constructs default names corresponding
%   to the shapefile names as provided by HydroSHEDS. The basin shapefiles
%   need to be located in the data/system folder of the CLIMADA flood
%   module.
%   res: resolution of shape file - the HydroSHEDS basin outline files are
%   provided in 15 arcsecond (res==15; default) and 30 arcsecond (res==30)
%   resolution.
%   check_plots: whether a plot of the centroids should be generated (=1),
%   or not (=0; default)
%   force_recalc: if set to 1, basin IDs are recalculated even if they
%   already exist (default is 0)
% OUTPUTS:
%   centroids: centroids with an additional field 'basin_ID' denoting the
%   river basin a centroid has been assigned to
%
% MODIFICATION HISTORY:
% Melanie Bieli, melanie.bieli@bluewin.ch, 20150301, initial
% Melanie Bieli, melanie.bieli@bluewin.ch, 20150314, add basin identification based on continents
% Gilles Stassen, gillesstassen@hotmail.com, 20150315, added use of subdir
% Gilles Stassen, gillesstassen@hotmail.com, 20150327, automatic download functionality
%-

% set global variables
global climada_global

% check input arguments
if ~climada_init_vars; return; end
if ~exist('centroids','var') || isempty(centroids)
    climada_centroids_load
end
if ~exist('res','var') || isempty(res), res = 15;  end
if ~exist('basin_shapefile','var') || isempty(basin_shapefile)
    basin_shapefile=''; end
if ~exist('check_plots', 'var')|| isempty(check_plots),
    check_plots = 0; end

% PARAMETERS
%
% Define bounding boxes of the areas covered by the HydroSHEDS basin
% shapefiles (hardwired, since not expected to change)
% In a bounding box [a,b,c,d], a refers to west, b to east, c to south,
% and d to north
bounding_box.ca = [-119, -60, 5, -39];  % Central America
bounding_box.na = [-138, -52, 24, 61];  % North America
bounding_box.sa = [-93, -32, -56, 15];  % South America
bounding_box.eu = [-14, 70, 12, 62];    % Europe
bounding_box.as = [57, 180, -12, 61];   % Asia
bounding_box.af = [-19, 55, -35, 38];   % Africa
bounding_box.au = [112, 180, -56, -10]; % Australia
bounding_box.names = {
    'Central America'
    'North America'
    'South America'
    'Europe'
    'Asia'
    'Africa'
    'Australia'
    };
%
% fields contains the continent abbreviations that can be used to construct
% file names
fields = fieldnames(bounding_box);
%
% default data directory
module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];
%
% resolution string (a component of the default shapefile names)
res_string = sprintf('%ds',res);
%-

% Define centroids bounding boxes
centroids_rect = [min(centroids.lon), max(centroids.lon),...
    min(centroids.lat) max(centroids.lat)];
centroids_edges_x = [centroids_rect(1),centroids_rect(2),...
    centroids_rect(2),centroids_rect(1)];
centroids_edges_y = [centroids_rect(3),centroids_rect(3),...
    centroids_rect(4),centroids_rect(4)];


% Load or read basin shapefile
% A) basin shapefile has been given as input
if ~isempty(basin_shapefile)
    exist_matfile = climada_check_matfile(basin_shapefile);
    if exist_matfile
        [fP,fN,~] = fileparts(basin_shapefile);
        basin_shapefile = [fP filesep fN '.mat'];
        load(basin_shapefile);
    elseif exist(basin_shapefile,'file')
        shapes = climada_shaperead(basin_shapefile);
    else
        cprintf([206 50 50]/255,'ERROR: did not find the specified file %s.\n',...
            basin_shapefile);
    end
    % B) use default basin shapefile
else
    % B.1): check whether centroids bounding box is contained in a
    % shapefile bounding box
    % Determine in which of the default HydroSHEDS shapefiles the
    % centroids are located
    fprintf('determining HydroSHEDS bounding box...\n')
    quit = 0;
    for i = 1:(length(fields)-1)
        box_edges_x = [bounding_box.(fields{i})(1), ...
            bounding_box.(fields{i})(2), bounding_box.(fields{i})(2),...
            bounding_box.(fields{i})(1)];
        box_edges_y = [bounding_box.(fields{i})(3), ...
            bounding_box.(fields{i})(3), bounding_box.(fields{i})(4),...
            bounding_box.(fields{i})(4)];
        in_box = inpolygon(centroids_edges_x, centroids_edges_y,...
            box_edges_x, box_edges_y);
        if isempty(in_box(in_box==0))
            quit=1;
            % found the centroids' corresponding basin shapefile
            shapefile_name = sprintf('%s_bas_%s_beta',fields{i},res_string);
            
            % look for file in root dir
            basin_dir = subdir([climada_global.root_dir filesep shapefile_name '.shp']);
            if ~isempty(basin_dir)
                basin_shapefile = basin_dir.name;
                % check if matfile already exists before reading the basin
                % shapefile
                exist_matfile = climada_check_matfile(basin_shapefile);
            else
                exist_matfile = 0;
            end
            
            if exist_matfile
                [fP,fN] = fileparts(basin_shapefile);
                basin_matfile = [fP filesep fN '.mat'];
                load(basin_matfile);
            elseif exist(basin_shapefile,'file')
                shapes = climada_shaperead(basin_shapefile);
            else
                % download directly
                fprintf('downloading HydroSHEDS basin shape files...')
                basin_shape_URL = ['http://earlywarning.usgs.gov/hydrodata/sa_shapefiles_zip/'...
                    shapefile_name '.zip'];
                basin_dir = [module_data_dir filesep 'system' filesep shapefile_name];
                basin_shapefile = unzip(basin_shape_URL,basin_dir);
                shapes = climada_shaperead(basin_shapefile,1);
                fprintf(' done\n')
            end % exist matfile
        end % isempty(in_box(in_box==0))
        if quit==1, break; end % no need to check the remaining boxes
    end % loop over fields
end % if ~isempty(basin_shapefile)

if isempty(basin_shapefile) && isfield(centroids, 'admin0_ISO3')
    % B.2): Centroids bounding box is not fully contained in one of
    % the shapefile bounding boxes. We therefore try to determine the
    % correct shapefile from continent names
    %         fprintf(['Centroids bounding box is not fully contained in '...
    %             'one of the \nHydroSHEDS shapefile bounding boxes.\n' ...
    %             'Using continent names instead...\n']);
    if exist(climada_global.map_border_file, 'file')
        load(climada_global.map_border_file)
    else
        cprintf([206 50 50]/255,['ERROR: climada_global.map_border_file '...
            'not found.\nPlease download the '...
            '<a href="https://github.com/davidnbresch/climada">'...
            'CLIMADA core module</a> from Github.\n']);
        return;
    end
    % Find country in the map border shapefile and determine the
    % continent
    
    [~,~,country_ndx] = climada_country_name(centroids.admin0_ISO3); % err msg in climada_country_name
    
    continent = shapes(country_ndx).CONTINENT;
    continent_ndx = strcmp(continent, bounding_box.names);
    if isempty(country_ndx)
        cprintf([206 50 50]/255,['ERROR: could not determine '...
            'default HydroSHEDS shapefile for %s\n', centroids.admin0_ISO3]);
        return;
    end
    shapefile_name = sprintf('%s_bas_%s_beta',...
        fields{continent_ndx},res_string);
    % special cases
    % Mexico would otherwise be assinged to North America
    if strcmp(centroids.admin0_name,'Mexico')
        shapefile_name = sprintf('ca_bas_%s_beta',res_string);
    end
    basin_shapefile = [module_data_dir filesep 'system' filesep ...
        shapefile_name '.shp'];
    
    % check if matfile already exists before reading the basin
    % shapefile
    exist_matfile = climada_check_matfile(basin_shapefile);
    if exist_matfile
        [fP,fN,~] = fileparts(basin_shapefile);
        basin_matfile = [fP filesep fN '.mat'];
        load(basin_matfile);
    elseif exist(basin_shapefile,'file')
        shapes = climada_shaperead(basin_shapefile,1);
    else
        fprintf('downloading HydroSHEDS basin shape files...')
        basin_shape_URL = ['http://earlywarning.usgs.gov/hydrodata/sa_shapefiles_zip/'...
            shapefile_name '.zip'];
        basin_dir = [module_data_dir filesep 'system' filesep shapefile_name];
        basin_shapefile = unzip(basin_shape_URL,basin_dir);
        shapes = climada_shaperead(basin_shapefile);
        fprintf(' done\n')
    end % exist_matfile
end % if isempty(basin_shapefile)

% Preselect basins such that only basins that overlap with the
% centroids structure will be considered for assignment of basin IDs
preselect_basin_IDs = false(length(shapes),1);

for basin_i = 1: length(shapes)
    
    if any(inpolygon(shapes(basin_i).X, shapes(basin_i).Y, ...
            centroids_edges_x, centroids_edges_y))
        preselect_basin_IDs(basin_i) = true;
    end
end

shapes(~preselect_basin_IDs) = [];

% Prepare input for basin_identify
temp = struct2cell(shapes).'; lon_polygons = temp(:,3); clear temp;
lon_polygons = lon_polygons';
temp = struct2cell(shapes).'; lat_polygons = temp(:,4); clear temp;
lat_polygons = lat_polygons';
temp = struct2cell(shapes).'; basin_names = temp(:,5); clear temp;
basin_names = basin_names';

% Assign basin ID to the centroids
fprintf('assigning basin IDs to %i centroids... ',length(centroids.centroid_ID));
basin_IDs = basin_identify(centroids.lon,centroids.lat,lon_polygons,...
    lat_polygons,basin_names);
centroids.basin_ID = basin_IDs;
fprintf('done \n')

% If required, generate a plot of the centroids highlighting the centroids
% that have been assigned a basin ID, as well as the basin outlines
if check_plots
    % calculate figure parameters
    markersize = 5;
    show_colorbar = 0;
    number_of_basins = length(unique(centroids.basin_ID));
    cmap = jet(number_of_basins);
    title_string = sprintf('Centroids in %s, colored by basin ID',...
        centroids.admin0_name);
    scale  = max(centroids.lon) - min(centroids.lon);
    ax_lim = [min(centroids.lon)-scale/30     max(centroids.lon)+scale/30 ...
        max(min(centroids.lat),-60)-scale/30  min(max(centroids.lat),80)+scale/30];
    
    figure('Name','Basin IDs','Color',[1 1 1]);
    % plot climada world borders
    climada_plot_world_borders(0.5);
    xlabel('Longitude'); ylabel('Latitude')
    axis(ax_lim)
    axis equal
    axis(ax_lim)
    
    % plot the centroids color-coded according to their basin IDs
    h=plotclr(centroids.lon, centroids.lat, centroids.basin_ID, '.',...
        markersize,show_colorbar,[],[],cmap);
    caxis([0 size(cmap,1)])
    title(title_string)
    
    %     unique_IDs = unique(centroids.basin_ID);
    %     cbar_label = unique_IDs(2:end);  % skipped 0, since not a valid basin ID
    %     for i = 1:length(number_of_basins)
    %         cbar_label(i) = unique_IDs(i+1);
    %     end
    %
    %     set(h,'YTick',0.5:1:size(cmap,1)-0.5,'yticklabel',cbar_label,'fontsize',12)
    
end % if check_plots

