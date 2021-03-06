function [centroids, shapes] = centroids_basin_ID(centroids, res, check_plots)
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
% Gilles Stassen, 20150408, clean up, the file finding piece in particular
%-

% set global variables
global climada_global

% check input arguments
if ~climada_init_vars; return; end
if ~exist('centroids','var') || isempty(centroids)
    climada_centroids_load
end
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
centroids_rect = [min(centroids.lon), max(centroids.lon),...
    min(centroids.lat) max(centroids.lat)];
centroids_edges_x = [centroids_rect(1),centroids_rect(2),...
    centroids_rect(2),centroids_rect(1)];
centroids_edges_y = [centroids_rect(3),centroids_rect(3),...
    centroids_rect(4),centroids_rect(4)];


% Load and read basin shapefile
shapefile_name = []; basin_shapefile = [];
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
    shapefile_name = sprintf('%s_bas_%s_beta',fields{continent_ndx},res_string);
    % special cases !!!!THERE ARE PROBABLY MORE THAN THE TWO BELOW!!!!
    % Mexico would otherwise be assinged to North America
    if strcmp(centroids.admin0_ISO3,'MEX') || strcmp(centroids.admin0_ISO3,'SLV')
        shapefile_name = sprintf('ca_bas_%s_beta',res_string);
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
            shapefile_name = sprintf('%s_bas_%s_beta',fields{i},res_string);
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

basin_shapefile = [module_data_dir filesep 'system' filesep shapefile_name filesep shapefile_name '.shp'];

[fP,fN] = fileparts(basin_shapefile);
basin_matfile = [fP filesep fN '.mat'];
if climada_check_matfile(basin_shapefile)
    fprintf('loading HydroSHEDS basin .mat file from %s \n',basin_matfile)
    load(basin_matfile);
elseif exist(basin_shapefile,'file')
    %fprintf('reading shape files from %s \n',basin_shapefile)
    shapes = climada_shaperead(basin_shapefile);
    %save(basin_matfile,'shapes')
    %fprintf('saved as .mat file to %s \n',basin_matfile)
else
    % download directly
    basin_shape_URL = ['http://earlywarning.usgs.gov/hydrodata/sa_shapefiles_zip/' shapefile_name '.zip'];
    fprintf('please wait - downloading HydroSHEDS basin shape files from %s \n',basin_shape_URL)
    basin_dir = [module_data_dir filesep 'system' filesep shapefile_name];
    basin_files = unzip(basin_shape_URL,basin_dir);
    for file_i=1:length(basin_files)
        [~,~,fE]=fileparts(basin_files{file_i});
        if strcmp(fE,'.shp')
            basin_shapefile = basin_files{file_i};
            break;
        end
    end
    shapes = climada_shaperead(basin_shapefile,1);
    %save(basin_matfile,'shapes')
    %fprintf('saved as .mat file to %s \n',basin_matfile)
end % exist matfile

% Preselect basins such that only basins that overlap with the
% centroids structure will be considered for assignment of basin IDs
preselect_basin_IDs = false(length(shapes),1);

for basin_i = 1: length(shapes)
    if any(inpolygon(shapes(basin_i).X, shapes(basin_i).Y, centroids_edges_x, centroids_edges_y))
        preselect_basin_IDs(basin_i) = true;
    end
end

shapes(~preselect_basin_IDs) = [];  % relevant shapes, also for arg out

if isempty(shapes) % incorrect shapefile found (probably)
    cprintf([1 0 0], 'ERROR: could not determine basin shapes for the centroids given\n')
    return
end

% Prepare input for basin_identify
tmp = struct2cell(shapes).'; lon_polygons = tmp(:,3); clear tmp;
lon_polygons = lon_polygons';
tmp = struct2cell(shapes).'; lat_polygons = tmp(:,4); clear tmp;
lat_polygons = lat_polygons';
tmp = struct2cell(shapes).'; basin_names = tmp(:,5); clear tmp;
basin_names = basin_names';

% Assign basin ID to the centroids
fprintf('assigning basin IDs to %i centroids... ',length(centroids.centroid_ID));
% loop over basins and check which centroids belong to them
for basin_i=1:length(basin_names)
   
    basin_str = basin_names{basin_i};
    tmp_lon = lon_polygons{basin_i};
    tmp_lat = lat_polygons{basin_i};
    tmp_lon = [tmp_lon tmp_lon(1)];   %to wrap around
    tmp_lat = [tmp_lat tmp_lat(1)];
        
    ndx = inpolygon(centroids.lon,centroids.lat,tmp_lon,tmp_lat);
    if(sum(ndx)>0)
        basin_IDs(ndx) = {basin_str};
    end
    clear ndx
    
end

basin_IDs(cellfun(@isempty,basin_IDs)) = {0};
basin_IDs = cell2mat(basin_IDs);

% basin_IDs = basin_identify(centroids.lon,centroids.lat,lon_polygons,lat_polygons,basin_names);
centroids.basin_ID = round(basin_IDs./100).*100;
fprintf('done \n')


% If required, generate a plot of the centroids highlighting the centroids
% that have been assigned a basin ID, as well as the basin outlines
if check_plots
    % calculate figure parameters
    markersize = 5;
    show_colorbar = 0;
    cmap = jet(length(unique(centroids.basin_ID)));
    title_string = sprintf('Centroids in %s, colored by basin ID',centroids.admin0_ISO3);
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

