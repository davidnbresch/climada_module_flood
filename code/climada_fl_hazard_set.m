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
if ~exist('check_plots',        'var'),     check_plots     = [];   end
if ~exist('fl_hazard_save_file',    'var'), fl_hazard_save_file = [];   end

% prompt for TR hazard event set if not given
if isempty(hazard_tr) % local GUI
    hazard_tr_file=[climada_global.data_dir filesep 'hazards' filesep '*.mat'];
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
    fl_hazard_save_file=[climada_global.data_dir filesep 'hazards' ...
        filesep 'TS_hazard.mat'];
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


basin_IDs                   =   unique(centroids.basin_ID);
number_of_basins            =   length(basin_IDs);

for basin_i = 1:number_of_basins
    c_ndx                       =   centroids.basin_ID == basin_IDs(basin_i);
    
    basin(basin_i).intensity    =   hazard_tr.intensity(:,c_ndx);
    basin(basin_i).lon          =   hazard_tr.lon(c_ndx); % would be the same as centroids.lon(c_ndx)
    basin(basin_i).lat          =   hazard_tr.lat(c_ndx);
    
    fl_score_sum                =   sum(centroids.flood_score(c_ndx));
    wet_index_sum               =   sum(centroids.wetness_index(c_ndx));
    
    if fl_score_sum ~= 0
        for event_i = 1 : hazard.event_count
            rain_sum                            =   sum(basin(basin_i).intensity(event_i,:),2);
            %hazard.intensity(event_i,c_ndx)     =   rain_sum .* (centroids.flood_score(c_ndx) ./ fl_score_sum);
            hazard.intensity(event_i,c_ndx)     =   rain_sum .* (centroids.wetness_index(c_ndx) ./ wet_index_sum);
        end 
    end
end
    

