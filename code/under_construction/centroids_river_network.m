function centroids = centroids_river_network(centroids,check_plots)
% identify river network in centroids
% MODULE:
%   flood
% NAME:
%	centroids_river_network
% PURPOSE:
%   Add rivers to centroids based on HydroSHEDS river network outline
%   shapefiles For more information and a technical documentation see
%   http://hydrosheds.org/page/hydrobasins,
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
%         .admin0_NAME    Country name
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
%   Gilles Stassen, gillesstassen@hotmail.com, 20150318
%-

% set global variables
global climada_global

% check input arguments
if ~climada_init_vars; return; end
if ~exist('centroids',   'var')|| isempty(centroids),climada_centroids_load; end
if ~exist('check_plots', 'var')|| isempty(check_plots),     check_plots = 0; end

% default data directory
module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];

if ~isfield(centroids,'admin0_NAME') || ~isfield(centroids,'admin0_ISO3') ...
        || ~isfield(centroids,'country_name')
    load(climada_global.map_border_file)
%     check_pts.X = [min(centroids.lon) min(centroids.lon) max(centroids.lon) max(centroids.lon) mean(centroids.lon)];
%     check_pts.Y = [min(centroids.lat) max(centroids.lat) min(centroids.lat) max(centroids.lat) mean(centroids.lat)];
    
    check_pts.X = [mean(centroids.lon)];
    check_pts.Y = [mean(centroids.lat)];

    %indexing
    shape_ndx = false(size(shapes));
    
    for shape_i = 1:length(shapes)
        % use rectangular shape for initial speed up
        bb_poly.X = [min(shapes(shape_i).X) min(shapes(shape_i).X) max(shapes(shape_i).X) max(shapes(shape_i).X)];
        bb_poly.Y = [min(shapes(shape_i).Y) max(shapes(shape_i).Y) max(shapes(shape_i).Y) min(shapes(shape_i).Y)];
        if any(inpolygon(check_pts.X,check_pts.Y,bb_poly.X,bb_poly.Y))
            shape_ndx(shape_i) = true;
        end
    end

%     admin0_ISO3 = shapes(shape_ndx).ADM0_A3;
    shapes = shapes(shape_ndx);
    % reset indexing
    shape_ndx = false(size(shapes));
    for shape_i = 1:length(shapes)
        % now use entire polygon
        if any(inpolygon(check_pts.X,check_pts.Y,shapes(shape_i).X,shapes(shape_i).Y))
            shape_ndx(shape_i) = true;
        end
    end
    shapes = shapes(shape_ndx);
    ISO3 = shapes.ADM0_A3;
else
    if isfield(centroids,'admin0_ISO3')
        ISO3 = centroids.admin0_ISO3;
    elseif isfield(centroids,'admin0_NAME')
        [~,ISO3,~] = climada_country_name(centroids.admin0_NAME);
    elseif isfield(centroids,'country_name')
        [~,ISO3,~] = climada_country_name(mode(centroids.country_name));
    end
end
if isempty(ISO3)
    cprintf([1 0 0 ], 'ERROR: unable to determine country ISO3 for given centroids \n')
end

% get relevant files from the web
rivers_shapefile = [module_data_dir filesep 'system' filesep ISO3 '_wat' filesep ISO3 '_water_lines_dcw.shp'];
rivers_dir = fileparts(rivers_shapefile);
if ~exist(rivers_shapefile,'file')
    rivers_URL = ['http://biogeo.ucdavis.edu/data/diva/wat/' ISO3 '_wat.zip'];
    fprintf('downloading rivers shapefiles from http://www.diva-gis.org/gdata ...');
    unzip(rivers_URL,rivers_dir);
    fprintf('done \n')
end        


shapes = climada_shaperead(rivers_shapefile,1);

fprintf('assigning river IDs to centroids...')
%init
centroids.river_ID = zeros(size(centroids.centroid_ID));
rivers.X = [];
rivers.Y = [];
rivers.ID =[];
for river_i = 1: length(shapes)
    for node_i = 1:length(shapes(river_i).X)
        dist_m = climada_geo_distance(shapes(river_i).X(node_i),shapes(river_i).Y(node_i),...
            centroids.lon,centroids.lat);
        [min_dist,min_dist_ndx] = min(dist_m);
        centroids.river_ID(min_dist_ndx) = river_i;
        if isnan(min_dist)
            rivers.X = [rivers.X NaN];
            rivers.Y = [rivers.Y NaN];
            rivers.ID =[rivers.ID river_i];
        else
            % for nicer plotting
            rivers.X = [rivers.X    centroids.lon(min_dist_ndx)];
            rivers.Y = [rivers.Y    centroids.lat(min_dist_ndx)];
            rivers.ID =[rivers.ID   river_i];
        end
    end
    rivers.X = [rivers.X    NaN];
    rivers.Y = [rivers.Y    NaN];
    rivers.ID =[rivers.ID   NaN];
end
fprintf(' done\n')

% just in case
centroids.river_shapes = shapes;

if check_plots
    climada_plot_world_borders;
    axis([min(centroids.lon) max(centroids.lon) min(centroids.lat) max(centroids.lat)])
    hold on 
%     for river_i = unique(centroids.river_ID(centroids.river_ID~=0)),
%         plot(centroids.lon(centroids.river_ID == river_i),...
%             centroids.lat(centroids.river_ID == river_i),'.b'); 
%     end
    plot(rivers.X,rivers.Y,'b')
end

