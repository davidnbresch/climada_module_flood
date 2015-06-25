function hazard = climada_fl_hazard_set(hazard_rf,centroids,hazard_set_file, check_plots)
% Generate flood hazard set from tr hazard set
% MODULE:
%   flood
% NAME:
%   climada_fl_hazard_set
% PURPOSE:
%   Generate flood hazard set from tr hazard set by distributing rainfall
%   volume according to flood scores (or wetness indices) of the centroids
% PREVIOUS STEP:
%   centroids_fl_prepare
% CALLING SEQUENCE:
%   hazard = climada_fl_hazard_set(hazard_tr,centroids,hazard_set_file, check_plots)
% EXAMPLE:
%   hazard = climada_fl_hazard_set('',centroids,hazard_set_file)
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
%   hazard_tr: a TR hazard event set, or a filename of a saved one
%       > prompted for if not given
%   hazard_set_file: the name of the hazard set file
%       > prompted for if not given
%   centroids: Climada centroids struct; the following fields are required:
%         .lat:           Latitude
%         .lon:           Longitude
%         .centroid_ID:   centroid ID
%         .admin0_name    Country name
% OUTPUTS:
%   hazard: a struct, the hazard event set, more for tests, since the
%       hazard set is stored as hazard_set_file, see code
%       lon(centroid_i): the longitude of each centroid
%       lat(centroid_i): the latitude of each centroid
%       centroid_ID(centroid_i): a unique ID for each centroid
%       peril_ID: just an ID identifying the peril
%       comment: a free comment, normally containing the time the hazard
%           event set has been generated
%       orig_years: the original years the set is based upon
%       orig_event_count: the original events
%       event_count: the total number of events in the set, including all
%           probabilistic ones, hence event_count>=orig_event_count
%       orig_event_flag(event_i): a flag for each event, whether it's an original
%           (1) or probabilistic (0) one
%       event_ID: a unique ID for each event
%       date: the creation date of the set
%       arr(event_i,centroid_i),sparse: the hazard intensity of event_i at
%           centroid_i
%       frequency(event_i): the frequency of each event
%       matrix_density: the density of the sparse array hazard.intensity
%       filename: the filename of the hazard event set (if passed as a
%           struct, this is often useful)
% MODIFICATION HISTORY:
% Gilles Stassen, gillesstassen@hotmail.com, 20150310, initial
% Melanie Bieli, melanie.bieli@bluewin.ch, 20150312, adjustments in rainfall distribution

hazard = []; % init

global climada_global
if ~climada_init_vars, return; end

% check arguments
if ~exist('hazard_rf',          'var'),     hazard_rf       = [];   end
if ~exist('centroids',          'var'),     centroids       = [];   end
if ~exist('check_plots',        'var'),     check_plots     = 0;    end
if ~exist('hazard_set_file',    'var'),     hazard_set_file = [];  	end

module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];

% prompt for TR hazard event set if not given
if isempty(hazard_rf) % local GUI
    hazard_tr_file=[climada_global.data_dir filesep 'hazards' filesep '*.mat'];
    [filename, pathname] = uigetfile(hazard_tr_file,...
        'Select a rainfall (eg TR) hazard event set:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        hazard_tr_file =fullfile(pathname,filename);
    end
    load(hazard_tr_file);
    hazard_rf=hazard;hazard=[];
end

% prompt for fl_hazard_save_file if not given
if isempty(hazard_set_file) % local GUI
    hazard_set_file=[climada_global.data_dir filesep 'hazards' filesep 'FL_hazard.mat'];
    [filename, pathname] = uiputfile(hazard_set_file, ...
        'Save new flood (FL) hazard event set as:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        hazard_set_file=fullfile(pathname,filename);
    end
end

hazard                  =   hazard_rf;
hazard.peril_ID         =   'FL';
hazard.date             =   datestr(now);
hazard.filename         =   hazard_set_file;
hazard.matrix_density   =   [];
hazard.comment          =   [];
hazard.intensity        =   zeros(size(hazard_rf.intensity));


if isfield(hazard,'rainfield_comment')
    hazard = rmfield(hazard, 'rainfield_comment');
end

if isempty(centroids)
    centroids.lon=hazard.lon';
    centroids.lat=hazard.lat';
    centroids.centroid_ID=hazard.centroid_ID;
end

centroids_n_flds = length(fieldnames(centroids));

centroids = climada_fl_centroids_prepare(centroids,0,0,'NO_SAVE');

% auto save if updated (i.e. new fields added)
if isfield(centroids,'filename') && length(fieldnames(centroids)) > centroids_n_flds
    fprintf('autosaving centroids with additional fields to %s \n',centroids.filename)
    save(centroids.filename,'centroids')
end
clear centroids_n_flds

if check_plots
    [x, y] = meshgrid(unique(centroids.lon),unique(centroids.lat));
    elev_tmp    = centroids.elevation_m; % init
%     elev_tmp(~centroids.onLand) = 0;
    gridded_elev = griddata(centroids.lon,centroids.lat,elev_tmp,x,y);
    figure('name',sprintf('Flood height in m'),'color','w')
    h = surf(x,y,gridded_elev, 'edgecolor','none');
    colormap(flipud(bone));
    hold on
end

if ~isfield(centroids,'basin_ID')
    centroids       = centroids_basin_ID(centroids);
end
if ~isfield(centroids,'elevation_m') && ~isempty(which('climada_read_srtm_DEM'))
    [~,centroids]   = climada_read_srtm_DEM('DL',centroids);
end   
if ~isfield(centroids,'TWI')
    centroids       = centroids_TWI(centroids);
end

basin_IDs           =   unique(centroids.basin_ID);
n_basins            =   length(basin_IDs);

if ~isfield(centroids,'onLand')
    if exist(climada_global.coastline_file,'file')
        fprintf('determining on land centroids... ')
        load(climada_global.coastline_file)
        onLand_ndx = inpolygon(centroids.lon,centroids.lat,shapes.X,shapes.Y);
        centroids.onLand( onLand_ndx)   = 1;
        centroids.onLand(~onLand_ndx)   = 0;
        fprintf('done \n')
    elseif exist(climada_global.map_border_file,'file')
        fprintf('determining on land centroids... ')
        load(climada_global.map_border_file)
        onLand_ndx = inpolygon(centroids.lon,centroids.lat,shapes.X,shapes.Y);
        centroids.onLand( onLand_ndx)   = 1;
        centroids.onLand(~onLand_ndx)   = 0;
        fprintf('done \n')
    else
        centroids.onLand                = 1;
        centroids.onLand(centroids.elevation_m < 0) = 0;
    end
end

% find minimum interval size
elev_res            =   min(diff(unique(centroids.elevation_m)));

% for progress mgmt

format_str	= '%s';
t0            = clock;

for basin_i = 1:n_basins
    % progress mgmt
    t_elapsed       = etime(clock,t0)/basin_i;
    n_remaining     = n_basins-basin_i;
    t_projected_sec = t_elapsed*n_remaining;
    if t_projected_sec<60
        msgstr = sprintf('est. %3.0f sec left (%i/%i basins): ',t_projected_sec, basin_i, n_basins);
    else
        msgstr = sprintf('est. %3.1f min left (%i/%i basins): ',t_projected_sec/60, basin_i, n_basins);
    end
    fprintf(format_str,msgstr);
    format_str = [repmat('\b',1,length(msgstr)) '%s'];
    
    % index of centroids belonging to basin_i
    c_ndx       =   (centroids.basin_ID == basin_IDs(basin_i)) & (centroids.onLand==1);
    
    % fl_score_sum                =   sum(centroids.flood_score(c_ndx));
    wet_index_sum               =   sum(centroids.TWI(c_ndx));
    
    n_events = hazard_rf.event_count;
    for event_i = 1 : n_events

        % index of rained-on centroids
        r_ndx   = hazard_rf.intensity(event_i,:) > 0;

        if ~any(r_ndx & c_ndx)
            continue
        end

        % index of floodable centroids, i.e. centroids lower in
        % elevation than highest r_ndx centroid
        fl_ndx  = c_ndx & (centroids.elevation_m <= max(centroids.elevation_m(r_ndx & c_ndx)));

        rain_sum  =   sum(hazard_rf.intensity(event_i,r_ndx & c_ndx),2);

        if wet_index_sum ~=0
                hazard.intensity(event_i,fl_ndx) = rain_sum .*(centroids.TWI(fl_ndx) ./ wet_index_sum);

        else
            hazard.intensity(event_i,fl_ndx)     =   rain_sum / sum(fl_ndx);
        end
    end
end
    


% flood only makes sense on land
% centroids field 'onLand':
%   1 (or higher, if more than one country):    on land,
%   0:                                          on sea,
%   and maximum value stands for buffer zone
%hazard.intensity(centroids.onLand==0) = 0;
%hazard.intensity(centroids.onLand==max(centroids.onLand))=0;
if hazard_set_file ~= 0
    fprintf('saving FL hazard set as %s\n',hazard_set_file);
    save(hazard_set_file,'hazard')
end

return