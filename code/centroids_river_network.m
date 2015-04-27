function centroids = centroids_river_network(centroids,res,check_plots)
% identify river network in centroids
% MODULE:
%   flood
% NAME:
%	centroids_river_network
% PURPOSE:
%   Add rivers to centroids based on HydroSHEDS river network outline
%   shapefiles For more information and a technical documentation see
%   http://hydrosheds.org/page/hydrobasins,
%   and to download the river network shapefiles see
%   http://hydrosheds.org/user/signin (requires registration)
%
% CALLING SEQUENCE:
%   centroids = centroids_river_network(river_shapes,centroids, check_plots)
% EXAMPLE:
%   centroids = centroids_river_network('',centroids,0)
% INPUTS:
%   centroids: Climada centroids struct; the following fields are required:
%         .lat:           Latitude
%         .lon:           Longitude
%         .centroid_ID:   centroid ID
%         .admin0_name    Country name
% OPTIONAL INPUT PARAMETERS:
%   river_shapes:   name (including full path) of the river network shapefile
%                   to be used. If set for 'DL', climada will
%                   automatically download the relevant files given the
%                   spatial extent of the centroids.
%   res:            resolution of shape file - the HydroSHEDS files are
%                   provided in 15 arcsecond (res==15; default) and 30
%                   arcsecond (res==30) resolution.
%   check_plots:    whether a plot of the centroids should be generated (=1),
%                   or not (=0; default)
% OUTPUTS:
%   centroids:      centroids with an additional field 'river' defining the
%                   river network as logicals.
% MODIFICATION HISTORY:
%   Gilles Stassen, gillesstassen@@hotmail.com, 20150318
%-

% set global variables
global climada_global

% check input arguments
if ~climada_init_vars; return; end
if ~exist('centroids',   'var')|| isempty(centroids),climada_centroids_load; end
if ~exist('res',         'var')|| isempty(res),             res = 15;        end
if ~exist('check_plots', 'var')|| isempty(check_plots),     check_plots = 0; end

if ~exist('res','var') || isempty(res), 
    res = 15;  
elseif ~ismember(res,[15 30])
    msg = sprintf('specified %i arcsecond resolution not available, choose 15 or 30 arcsecond resolution: ',res);
    res = input(msg);
    if ~ismember(res,[15 30])
        cprintf([1 0 0],'Oohh come on! What did I tell you? - ABORTING\n')
        return
    end
end
if ~exist('check_plots', 'var')|| isempty(check_plots),check_plots = 0; end

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
centroids_rect = [min(centroids.lon), max(centroids.lon),min(centroids.lat) max(centroids.lat)];
centroids_edges_x = [centroids_rect(1),centroids_rect(2),centroids_rect(2),centroids_rect(1)];
centroids_edges_y = [centroids_rect(3),centroids_rect(3),centroids_rect(4),centroids_rect(4)];


% Load and read basin shapefile
shapefile_name = []; river_shapefile = [];
if isfield(centroids, 'admin0_ISO3')
    % use continent name to select correct basin shape file
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
    shapefile_name = sprintf('%s_riv_%s_beta',fields{continent_ndx},res_string);
    % special cases
    % Mexico would otherwise be assinged to North America
    if strcmp(centroids.admin0_ISO3,'MEX') || strcmp(centroids.admin0_ISO3,'SLV')
        shapefile_name = sprintf('ca_riv_%s',res_string);
    end 
else
    % find shapefile that encloses centroids
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
            % found the centroids' corresponding basin shapefile
            shapefile_name = sprintf('%s_riv_%s',fields{i},res_string);
            break
        end 
    end
end

if isempty(shapefile_name) % shapefiles not found
    cprintf([1 0 0], 'ERROR: could not determine basin shapes for the centroids given\n')
    return
end

% look for file in root dir
% basin_dir = subdir([fileparts(climada_global.root_dir) filesep shapefile_name '.shp']);
% if ~isempty(basin_dir)
%     basin_shapefile = basin_dir.name;
% end

river_shapefile = [module_data_dir filesep 'system' filesep shapefile_name filesep shapefile_name '.shp'];

[fP,fN] = fileparts(river_shapefile);
river_matfile = [fP filesep fN '.mat'];
if climada_check_matfile(river_shapefile)
    fprintf('loading .mat file from %s \n',river_matfile)
    load(river_matfile);
elseif exist(river_shapefile,'file')
    %fprintf('reading shape files from %s \n',basin_shapefile)
    shapes = climada_shaperead(river_shapefile);
    %save(basin_matfile,'shapes')
    %fprintf('saved as .mat file to %s \n',basin_matfile)
else
    % download directly
    river_shape_URL = ['http://earlywarning.usgs.gov/hydrodata/sa_shapefiles_zip/' shapefile_name '.zip'];
    fprintf('please wait - downloading HydroSHEDS basin shape files from %s \n',river_shape_URL)
    river_dir = [module_data_dir filesep 'system' filesep shapefile_name];
    river_files = unzip(river_shape_URL,river_dir);
    for file_i=1:length(river_files)
        [~,~,fE]=fileparts(river_files{file_i});
        if strcmp(fE,'.shp')
            river_shapefile = river_files{file_i};
            break;
        end
    end
    shapes = climada_shaperead(river_shapefile,1);
    %save(basin_matfile,'shapes')
    %fprintf('saved as .mat file to %s \n',basin_matfile)
end % exist matfile

% Preselect basins such that only basins that overlap with the
% centroids structure will be considered for assignment of basin IDs
preselect_rivers = false(length(shapes),1);

if check_plots
    climada_plot_world_borders;
    axis([min(centroids.lon) max(centroids.lon) min(centroids.lat) max(centroids.lat)])
    hold on
end

for river_i = 1: length(shapes)
    if any(inpolygon(shapes(river_i).X, shapes(river_i).Y, centroids_edges_x, centroids_edges_y))
        preselect_rivers(river_i) = true;
        if check_plots
            plot(shapes(river_i).X,shapes(river_i).Y,'b')
        end
    end
end
hold off

shapes(~preselect_rivers) = [];



% If required, generate a plot of the centroids highlighting the centroids
% that have been assigned a basin ID, as well as the basin outlines
if check_plots
    % calculate figure parameters
    markersize = 5;
    show_colorbar = 0;
    number_of_basins = length(unique(centroids.basin_ID));
    cmap = jet(number_of_basins);
    title_string = sprintf('Centroids in %s, colored by basin ID',centroids.country_name{1});
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

