function hazard = climada_ls_hazard_trigger(centroids,n_events,...
    wiggle_factor_TWI,condition_TWI, wiggle_factor_slope,condition_slope)

% Assessing of the susceptibility of shallow landslides. Triggering areas of landslides (1/0) are 
% for e.g. 100 events.
% MODULE:
%   flood
% NAME:
%   climada_ls_hazard_trigger
% PURPOSE:
%   
% CALLING SEQUENCE:
%   
% EXAMPLE:
%   
% INPUTS:
%   
% OPTIONAL INPUT PARAMETERS:
%  
% OUTPUTS:
%  
% MODIFICATION HISTORY:
% Thomas Rölli, thomasroelli@gmail.com, 20180201, init

global climada_global
if ~climada_init_vars, return; end

% check arguments
if ~exist('centroids', 'var'), centroids = []; end
if ~exist('n_events', 'var'), n_events = []; end
if ~exist('hazard_set_file', 'var'), hazard_set_file = []; end
if ~exist('wiggle_factor', 'var'), wiggle_factor_TWI = []; end
if ~exist('TWI_condition', 'var'), condition_TWI = []; end
if ~exist('wiggle_factor_slope', 'var'), wiggle_factor_slope = []; end
if ~exist('slope_condition', 'var'), condition_slope = []; end

% PARAMETERS
if isempty(n_events); n_events = 100; end
if isempty(wiggle_factor_TWI); wiggle_factor_TWI = 0.35; end
if isempty(condition_TWI); condition_TWI = 0.95; end
if isempty(wiggle_factor_slope); wiggle_factor_slope = 0.2; end
if isempty(condition_slope); condition_slope = 0.45; end

%initiate hazard 
hazard = []; % init
hazard.lon = centroids.lon;
hazard.lat = centroids.lat;
hazard.centroid_ID = 1:numel(hazard.lon);
hazard.peril_ID         = 'LS';
hazard.orig_years       = 10000;
hazard.orig_event_count = n_events;
hazard.event_count      = n_events;
hazard.event_ID         = 1:n_events;
hazard.orig_event_flag  = ones(1,n_events);
hazard.yyyy             = ones(1,n_events);
hazard.mm               = ones(1,n_events);
hazard.dd               = ones(1,n_events);
hazard.intensity        = sparse(n_events,numel(hazard.lon));
hazard.name             = cell(1,n_events);
hazard.frequency        = ones(1,n_events)/hazard.orig_years;
hazard.comment          = centroids.comment;
hazard.date             = datestr(now);
hazard.units            = 'binary';
hazard.orig_yearset     = [];
hazard.filename         = hazard_set_file;
hazard.matrix_density   =  0.01; % estimate

% create TWI wiggle as the sum of TWI_norm + TWI_delta
% create slope wiggle as the sum of slope_factor + slope_delta
% TWI_wiggle = zeros(size(hazard.intensity)); %init
% slope_wiggle = zeros(size(hazard.intensity)); %init

delta_TWI = rand(size(hazard.intensity)) * wiggle_factor_TWI;
delta_slope = rand(size(hazard.intensity)) * wiggle_factor_slope;

for e_i = 1:n_events
    wiggle_TWI = centroids.TWI_norm + delta_TWI(e_i,:);
    wiggle_slope = centroids.slope_factor + delta_slope(e_i,:);
        
    % check where landslides occur
    ls_occurence(e_i,:) = wiggle_TWI>condition_TWI & wiggle_slope>condition_slope ;
    
end

% create sparse matrix
hazard.intensity = sparse(ls_occurence);
hazard.fraction = spones(hazard.intensity);

end