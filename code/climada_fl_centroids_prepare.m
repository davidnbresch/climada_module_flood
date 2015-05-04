function centroids = climada_fl_centroids_prepare(centroids, check_plots, force_recalc, save_file)
% Equip centroids with all necessary hydrological information 
% MODULE:
%   flood
% NAME:
%	climada_fl_centroids_prepare
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
% Gilles Stassen, gillesstassen@hotmail.com, 20150414, cleanup
% Gilles Stassen, 20150416, added save functionality, added BD, SD
%-

global climada_global

% check input arguments
if ~climada_init_vars; return; end
if ~exist('centroids','var') || isempty(centroids)
    climada_centroids_load
end
if ~exist('check_plots', 'var') || isempty(check_plots)
    check_plots = 0;        end
if ~exist('force_recalc','var') || isempty(force_recalc), 
    force_recalc = 0;       end
if ~exist('save_file',   'var'),    save_file = 'AUTO';     end

% Step 1: Compute centroid elevation fro
if ~isfield(centroids,'elevation_m')
    [~,centroids] = climada_read_srtm_DEM('DL',centroids,[],[],check_plots);
end

% Step 2: Calculate topographic wetness index (TWI) 
if ~isfield(centroids,'TWI') || force_recalc
    centroids = centroids_TWI(centroids, check_plots);
end

% Step 3: Delineate basins
if ~isfield(centroids,'basin_ID') || force_recalc
    res = 15; % hardwired
    centroids = centroids_basin_ID(centroids, res, check_plots);
end

% Step 3: Compute daily evapotranspiration (ET)
if ~isfield(centroids,'ET_mm_day') || force_recalc
    centroids = centroids_ET(centroids, check_plots);
end

% % Step 4: Assign soil wetness index (SWI)
% if ~isfield(centroids,'SWI') || force_recalc
%     centroids = centroids_SWI(centroids, check_plots);
% end

% Step 5: Assign available water-holding capacity of the soil (WHC)
if ~isfield(centroids,'WHC_mm') || force_recalc
    centroids = centroids_WHC(centroids, check_plots);
end

% Step 6: Assign soil bulk density values (BD)
if ~isfield(centroids, 'BD_kg_m3')
    centroids = centroids_BD(centroids,check_plots);
end

% Step 7: add soil depth
if ~isfield(centroids, 'SD_m')
    centroids = centroids_SD(centroids,check_plots);
end

% Step 8: add leaf area index
if ~isfield(centroids, 'LAI')
    centroids = centroids_LAI(centroids,check_plots);
end

% save the new centroids struct
if strcmp(save_file,'AUTO')
    if isfield(centroids, 'filename')
        fprintf('autosaving centroids as %s \n',centroids.filename)
        save(centroids.filename, 'centroids')
    else
        if isfield(centroids, 'admin0_ISO3')
            ISO3 = centroids.admin0_ISO3;
            save_file = [climada_global.data_dir filesep 'system' filesep 'centroids_' ISO3 '_' datestr(now,'ddmmyy') '.mat'];
        else
            save_file = [climada_global.data_dir filesep 'system' filesep 'centroids_' datestr(now,'ddmmyy') '.mat'];
        end
        fprintf('autosaving centroids as %s \n', save_file)
        centroids.filename = save_file;
        save(save_file, 'centroids')
    end
elseif ischar(save_file) && ~strcmp(save_file,'NO_SAVE')
    fprintf('saving centroids as %s \n', save_file)
    centroids.filename = save_file;
    save(save_file, 'centroids')        
end

