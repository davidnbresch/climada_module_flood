function centroids = fl_centroids_prepare(centroids, res, basin_shapefile, check_plots, force_recalc)
% Equip centroids with all necessary hydrological information 
% MODULE:
%   flood
% NAME:
%	fl_centroids_prepare
% PURPOSE:
%   Equip centroids with flood scores, wetness indices, basin IDs, 
%   evapotranspiration values, soil wetness indices, and values for the 
%   available water-holding capacity of the soil, such that the centroids
%   then contain all the information needed for the generation of 
%   a flood hazard set. 
%       In step 1, flood scores and topographic wetness indices are 
%       calculated for the given centroids 
%       (see function centroids_fl_score_calc for details)
%       In step 2, basin IDs are assigned to the centroids (see function
%       centroids_basinID_assign for details)
%       In step 3, evapotranspiration (ET) is calculated for the centroids
%       (see centroids_ET_assign for details)
%       In step 4, soil wetness index (SWI) is calculated for the centroids 
%       (see centroids_SWI_assign for details)
%       In step 5, available water-holding capacity of the soil (WHC) is 
%       calculated for the centroids (see centroids_WHC_assign for details)
% CALLING SEQUENCE:
%   centroids = fl_centroids_prepare(centroids, res, basin_shapefile,...
%       check_plots, force_recalc)
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
%   centroids: centroids with three additional fields: 
%       centroids.flood_score: flow accumulation 
%       centroids.topo_wetness_index: topographic wetness index 
%       centroids.basin_ID: river basins 
%       centroids.evapotranspiration: evapotranspiration (mm/yr)
%       centroids.soil_wetness_index: soil wetness index (%)
%       centroids.water_holding_capacity: available water-holding capacity
%           of the soil (mm)
%  
% MODIFICATION HISTORY:
% Melanie Bieli, melanie.bieli@bluewin.ch, 20150311, initial
% Melanie Bieli, melanie.bieli@bluewin.ch, 20150319, added LAI
% Melanie Bieli, melanie.bieli@bluewin.ch, 20150321, added ET
% Melanie Bieli, melanie.bieli@bluewin.ch, 20150324, added FC
% Melanie Bieli, melanie.bieli@bluewin.ch, 20150325, added SWI and WHC;
% removed variables that will presumably not be used 

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
    basin_shapefile='';     end
if ~exist('check_plots', 'var') || isempty(check_plots)
    check_plots = 0;        end
if ~exist('force_recalc','var') || isempty(force_recalc), 
    force_recalc = 0;       end

fprintf('Preparing %i centroids...\n', length(centroids.centroid_ID))
% Step 1: Calculate flood scores 
centroids = centroids_fl_score_calc(centroids, check_plots, force_recalc);

% Step 2: Assign basin IDs
centroids = centroids_basinID_assign(centroids, res, basin_shapefile, ...
    check_plots, force_recalc);

% Step 3: Assign evapotranspiration (ET)
centroids = centroids_ET_assign(centroids, check_plots, force_recalc);

% Step 4: Assign soil wetness index (SWI)
centroids = centroids_SWI_assign(centroids, check_plots, force_recalc);

% Step 5: Assign available water-holding capacity of the soil (WHC)
centroids = centroids_WHC_assign(centroids, check_plots, force_recalc);
end



