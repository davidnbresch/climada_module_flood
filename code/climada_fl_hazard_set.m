function hazard = climada_fl_hazard_set(hazard_tr,centroids,fl_hazard_save_file, check_plots)
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
if ~exist('hazard_tr',          'var'),     hazard_tr       = [];   end
if ~exist('centroids',          'var'),     centroids       = [];   end
if ~exist('check_plots',        'var'),     check_plots     = 0;   end
if ~exist('fl_hazard_save_file',    'var'), fl_hazard_save_file = [];   end

module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];

% prompt for TR hazard event set if not given
if isempty(hazard_tr) % local GUI
    hazard_tr_file=[module_data_dir filesep 'hazards' filesep '*.mat'];
    [filename, pathname] = uigetfile(hazard_tr_file,...
        'Select a rainfall hazard event set:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        hazard_tr_file =fullfile(pathname,filename);
    end
    load(hazard_tr_file);
end

% prompt for fl_hazard_save_file if not given
if isempty(fl_hazard_save_file) % local GUI
    fl_hazard_save_file=[module_data_dir filesep 'hazards' ...
        filesep 'test_FL_hazard.mat'];
    [filename, pathname] = uiputfile(fl_hazard_save_file, ...
        'Save new rainfall hazard event set as:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        fl_hazard_save_file=fullfile(pathname,filename);
    end
end

hazard                  =   hazard_tr;
hazard.peril_ID         =   'FL';
hazard.date             =   datestr(now);
hazard.filename         =   fl_hazard_save_file;
hazard.matrix_density   =   [];
hazard.comment          =   [];
hazard.intensity        =   zeros(size(hazard_tr.intensity));


if isfield(hazard,'rainfield_comment')
    hazard = rmfield(hazard, 'rainfield_comment');
end

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
    centroids = centroids_basinID_assign(centroids);
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

format_str_b	= '%s';
t0_b            = clock;

for basin_i = 1:n_basins
    % progress mgmt
    t_elapsed       = etime(clock,t0_b)/basin_i;
    n_remaining     = n_basins-basin_i;
    t_projected_sec = t_elapsed*n_remaining;
    msgstr          = ''; 
    if t_projected_sec<60
        msgstr_b = sprintf('est. %3.0f sec left (%i/%i basins): ',t_projected_sec, basin_i, n_basins);
    else
        msgstr_b = sprintf('est. %3.1f min left (%i/%i basins): ',t_projected_sec/60, basin_i, n_basins);
    end
    fprintf(format_str_b,msgstr_b);
    
    % index of centroids belonging to basin_i
    c_ndx       =   (centroids.basin_ID == basin_IDs(basin_i));% & (centroids.onLand==1);
    
    % fl_score_sum                =   sum(centroids.flood_score(c_ndx));
    wet_index_sum               =   sum(centroids.topo_wetness_index(c_ndx));
    
        % for progress mgmt
        mod_step    = 10;
        format_str	= '%s';
        t0          = clock;
        n_events    = hazard.event_count;
        for event_i = 1 : n_events
            
            % index of rained-on centroids
            r_ndx   = hazard_tr.intensity(event_i,:) > 0;
            
            if ~any(r_ndx & c_ndx)
                continue
            end
            
            % index of floodable centroids, i.e. centroids lower in
            % elevation than highest r_ndx centroid
            fl_ndx  = c_ndx & (centroids.elevation_m <= max(centroids.elevation_m(r_ndx & c_ndx)));
            
            rain_sum                            =   sum(hazard_tr.intensity(event_i,r_ndx & c_ndx),2);
            %hazard.intensity(event_i,fl_ndx)     =   rain_sum .* (centroids.flood_score(fl_ndx) ./ fl_score_sum);
            if wet_index_sum ~=0
                hazard.intensity(event_i,fl_ndx)     =   rain_sum .* (centroids.topo_wetness_index(fl_ndx) ./ wet_index_sum);
            else
                hazard.intensity(event_i,fl_ndx)     =   rain_sum / sum(fl_ndx);
            end
            
            % progress mgmt
            if mod(event_i,mod_step)==0
                mod_step        = 100;
                n_remaining     = n_events-event_i;
                t_projected_sec = t_elapsed*n_remaining;
                msgstr          = sprintf('%i/%i events', event_i, n_events);
                fprintf(format_str,msgstr);
                format_str=[repmat('\b',1,length(msgstr)) '%s'];
             end
        end
        format_str_b = [repmat('\b',1,length(msgstr_b)+length(msgstr)) '%s'];
end
    
return

% flood only makes sense on land
% centroids field 'onLand':
%   1 (or higher, if more than one country):    on land,
%   0:                                          on sea,
%   and maximum value stands for buffer zone
hazard.intensity(centroids.onLand==0) = 0;
hazard.intensity(centroids.onLand==max(centroids.onLand))=0;

fprintf('saving FL hazard set as %s\n',fl_hazard_save_file);
save(fl_hazard_save_file,'hazard')
