function centroids = centroids_basinID_assign(centroids, res, basin_shapefile, check_plots, force_recalc)
% assign basin ID to centroids
% MODULE:
%   tbd
% NAME:
%	country_basinID_assign
% PURPOSE:
%   Assign basin IDs to centroids based on HydroSHEDS basin outline 
%   shapefiles For more information and a technical documentation see 
%   http://hydrosheds.org/page/hydrobasins, 
%   and to download the basin shapefiles see 
%   http://hydrosheds.org/user/signin (requires registration) 
%   
% CALLING SEQUENCE:
%   centroids_basinID = centroids_basinID_assign(basin_shapefile,centroids, check_plots)
% EXAMPLE:
%   centroids_basinID = centroids_basinID_assign('',centroids,0)
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
% TO DO:
%   Add code to determine the basin shapefile covering the given centroids
%   from continent names (in case the bounding box approach does not work,
%   i.e. the centroids bounding box is not fully contained in one of the
%   shapefile bounding boxes)
% MODIFICATION HISTORY:
% Melanie Bieli, melanie.bieli@bluewin.ch, 20150225, initial


% set global variables
global climada_global

% check input arguments 
if ~climada_init_vars; return; end
if ~exist('centroids','var') || isempty(centroids)
    cprintf([206 50 50]/255,['Error: Missing input ''centroids'' to '...
        'function centroids_basin_ID_assign, can''t proceed.\n']);
    return;
end
if ~exist('res','var') || isempty(res), res = 15;  end
if ~exist('basin_shapefile','var') || isempty(basin_shapefile)
    basin_shapefile=''; end
if ~exist('check_plots', 'var')|| isempty(check_plots), 
    check_plots = 0; end
if ~exist('force_recalc','var') || isempty(force_recalc), 
    force_recalc = 0;  end


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

% We only calculate the basin IDs if the centroids do not come equipped
% with them
if ~isfield(centroids,'basin_ID') || force_recalc
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
            cprintf([206 50 50]/255,'Error: Did not find file %s.\n',...
                basin_shapefile);
        end
    % B) use default basin shapefile
    else
        % B.1): check whether centroids bounding box is contained in a
        % shapefile bounding box
        % Determine in which of the default HydroSHEDS shapefiles the 
        % centroids are located
        fprintf('Determining HydroSHEDS bounding box...\n')
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
                shapefile_name = sprintf('%s_bas_%s_beta.shp',...
                    fields{i},res_string);
                basin_shapefile = [module_data_dir filesep 'system' filesep ...
                    shapefile_name];
                % check if matfile already exists before reading the basin
                % shapefile
                exist_matfile = climada_check_matfile(basin_shapefile);
                if exist_matfile
                    [fP,fN,~] = fileparts(basin_shapefile);
                    basin_matfile = [fP filesep fN '.mat'];
                    load(basin_matfile);
                elseif exist(basin_shapefile,'file')
                    shapes = climada_shaperead(basin_shapefile);
                else
                    cprintf([206 50 50]/255,['Error: Missing basin shapefile %s.\n'...
                        'Please <a href="http://hydrosheds.org/user/signin">'...
                        'download it from the HydroSHEDS database</a>' ...
                        '(requires registration).\n'], shapefile_name);
                    return;
                end % exist matfile
            end % isempty(in_box(in_box==0))
            if quit==1, break; end % no need to check the remaining boxes
        end % loop over fields     
    end % if ~isempty(basin_shapefile)
    
    if isempty(basin_shapefile)
        % B.2): Centroids bounding box is not fully contained in one of
        % the shapefile bounding boxes. We therefore try to determine the
        % correct shapefile from continent names
        fprintf(['Centroids bounding box is not fully contained in '...
            'one of the \nHydroSHEDS shapefile bounding boxes.\n' ...
            'Using continent names instead...\n']);
        if exist(climada_global.map_border_file, 'file')
            load(climada_global.map_border_file)
        else
            cprintf([206 50 50]/255,['\nError: climada_global.map_border_file '...
                'not found.\nPlease download the '...
                '<a href="https://github.com/davidnbresch/climada">'...
                'CLIMADA core module</a> from Github.\n']);
            return;
        end
        % Find country in the map border shapefile and determine the
        % continent
        country_names = {shapes.NAME}.';
        if isfield(centroids,'admin0_name')
            country_index = find(strcmp(centroids.admin0_name,country_names));
        else
            cprintf([206 50 50]/255,['\nError: Could not determine'...
                'default HydroSHEDS shapefile for %s \n',centroids]);
            return;
        end
        if isempty(country_index)
            cprintf([206 50 50]/255,['\nError: Could not determine'...
                'default HydroSHEDS shapefile for %s\n', centroids.admin0_name]);
            return;
        end
        continent = shapes(country_index).CONTINENT;
        continent_index = find(strcmp(continent, bounding_box.names));
        if isempty(country_index)
            cprintf([206 50 50]/255,['\nError: Could not determine'...
                'default HydroSHEDS shapefile for %s\n', centroids.admin0_name]);
            return;
        end
        shapefile_name = sprintf('%s_bas_%s_beta.shp',...
            fields{continent_index},res_string);
        basin_shapefile = [module_data_dir filesep 'system' filesep ...
                    shapefile_name];
        % check if matfile already exists before reading the basin
        % shapefile
        exist_matfile = climada_check_matfile(basin_shapefile);
        if exist_matfile
            [fP,fN,~] = fileparts(basin_shapefile);
            basin_matfile = [fP filesep fN '.mat'];
            load(basin_matfile);
        elseif exist(basin_shapefile,'file')
            shapes = climada_shaperead(basin_shapefile);
        else
            cprintf([206 50 50]/255,['Error: Missing basin shapefile %s.\n'...
                'Please <a href="http://hydrosheds.org/user/signin">'...
                'download it from the HydroSHEDS database</a>' ...
                '(requires registration).\n'], shapefile_name);
            return;
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
    fprintf('Identifying basins for %i centroids...\n',...
        length(centroids.centroid_ID));
    basin_IDs = basin_identify(centroids.lon,centroids.lat,lon_polygons,...
        lat_polygons,basin_names);
    centroids.basin_ID = basin_IDs;
    cprintf([23 158 58]/255,['Successfully completed calculation of '...
        'basin IDs.\n']);
else
    % centroids already have a field 'basin_ID'
    cprintf([23 158 58]/255,'Skipped - centroids already have basin IDs.\n')
end % if ~isfield(centroids,'basin_ID')

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

end

