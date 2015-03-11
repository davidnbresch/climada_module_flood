function centroids = fl_centroids_prepare(centroids, res, basin_shapefile, check_plots, force_recalc)
% Assign flood score and basin IDs to centroids
% MODULE:
%   tbd
% NAME:
%	fl_centroids_prepare
% PURPOSE:
%   Equip centroids with flood scores and basin IDs, such that they can be
%   used for the generation of a flood hazard set. 
%       In step 1, flood scores are calculated for the given centroids (see
%       function centroids_fl_score_calc for details)
%       In step 2, basin IDs are assigned to the centroids (see function
%       centroids_basinID_assign for details)
% CALLING SEQUENCE:
%   centroids = fl_centroids_prepare(centroids, res, basin_shapefile, check_plots)
% EXAMPLE:
%   centroids = fl_centroids_prepare(centroids,15,'',1)
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
%   provided in 15 arcsecond (res==15) and 30 arcsecond (res==30; default)
%   resolution.
%   check_plots: whether a set of plots showing topography (elevation), 
%   slope, aspect angle, flood scores and basins should be generated (=1),
%   or not (=0; default)
%   force_recalc: if set to 1, recalculate flood scores and basin IDs, even
%   if they already exist (default is 0)
% OUTPUT:
%   centroids: centroids with two additional fields: 
%       centroids.flood_score: flow accumulation for each centroid
%       centroids.basin_ID: river basins the centroid have been assigned to
%  
% MODIFICATION HISTORY:
% Melanie Bieli, melanie.bieli@bluewin.ch, 20150311, initial

global climada_global

% check input arguments
if ~climada_init_vars; return; end
if ~exist('res','var') || isempty(res), res = 15;  end
if ~exist('centroids','var') || isempty(centroids)
    cprintf([206 50 50]/255,['Error: Missing input ''centroids'' to '...
        'function centroids_fl_score_calc, can''t proceed.\n']);
    return; 
end
if ~exist('basin_shapefile','var') || isempty(basin_shapefile)
    basin_shapefile=''; end
if ~exist('check_plots', 'var') || isempty(check_plots)
    check_plots = 0;    end
if ~exist('force_recalc','var') || isempty(force_recalc), 
    force_recalc = 0;   end

fprintf('Preparing %i centroids...\n', length(centroids.centroid_ID))
% Step 1: Calculate flood scores 
centroids = centroids_fl_score_calc(centroids, check_plots);

% Step 2: Assign basin IDs
centroids = centroids_basinID_assign(centroids, res, basin_shapefile, ...
    check_plots, force_recalc);
end



