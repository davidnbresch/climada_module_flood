function hazard = climada_fl_hazard_set(hazard_tr,centroids,hazard_set_file, check_plots)

hazard = []; % init

global climada_global
if ~climada_init_vars, return; end

if ~exist('hazard_tr',          'var'),     hazard_tr       = [];   end
if ~exist('centroids',          'var'),     centroids       = [];   end
if ~exist('check_plots',        'var'),     check_plots     = [];   end
if ~exist('hazard_set_file',    'var'),     hazard_set_file = [];   end
% prompt for TR hazard event set if not given
if isempty(hazard_tr) % local GUI
    tr_hazard_set_file=[climada_global.data_dir filesep 'hazards' filesep '*.mat'];
    [filename, pathname] = uigetfile(tr_hazard_set_file, 'Select a rainfall hazard event set:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        tr_hazard_set_file=fullfile(pathname,filename);
    end
    load(tr_hazard_set_file);
end

% % prompt for hazard_set_file if not given
% if isempty(hazard_set_file) % local GUI
%     hazard_set_file=[climada_global.data_dir filesep 'hazards' filesep 'TS_hazard.mat'];
%     [filename, pathname] = uiputfile(hazard_set_file, 'Save new rainfall hazard event set as:');
%     if isequal(filename,0) || isequal(pathname,0)
%         return; % cancel
%     else
%         hazard_set_file=fullfile(pathname,filename);
%     end
% end

hazard                  =   hazard_tr;
hazard.peril_ID         =   'FL';
hazard.date             =   datestr(now);
hazard.filename         =   hazard_set_file;
hazard                  =   rmfield(hazard, 'rainfield_comment');
hazard.matrix_density   =   [];
hazard.comment          =   [];
hazard.intensity        =   zeros(size(hazard_tr.intensity));


for basin_i = centroids.no_basins
    basin_IDs                   =   unique(centroids.basin_ID);
    c_ndx                       =   centroids.basin_ID == basin_IDs(basin_i);
    
    basin(basin_i).intensity    =   hazard_tr.intensity(:,c_ndx);
    basin(basin_i).lon          =   hazard_tr.lon(c_ndx);
    basin(basin_i).lat          =   hazard_tr.lat(c_ndx);
    
    fl_score_sum                =   sum(centroids.FL_score(c_ndx));
    
    if fl_score_sum ~= 0
        for event_i = 1 : hazard.event_count
            rain_sum                            =   sum(basin(basin_i).intensity(event_i,:),2);
            hazard.intensity(event_i,c_ndx)     =   rain_sum .* (centroids.FL_score(c_ndx) ./ fl_score_sum);
        end 
    end
end
    

