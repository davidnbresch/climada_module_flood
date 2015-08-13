function [hazard, centroidsORhazard] = climada_ls_hazard_set(hazardORcentroids, hazard_set_file, focus_areas, check_plots, chrono_check, bool_hazard_check)
% Generate ls hazard set from rainfall hazard set
% MODULE:
%   flood
% NAME:
%   climada_ls_hazard_set
% PURPOSE:
%   Generate land slide hazard set from tr hazard set by distributing rainfall
%   volume according to wetness indices of the centroids
% PREVIOUS STEP:
%   centroids_fl_prepare
% CALLING SEQUENCE:
%   hazard = climada_ls_hazard_set(hazard_rf,centroids,hazard_set_file, check_plots)
% EXAMPLE:
%   hazard = climada_ls_hazard_set('',centroids,hazard_set_file)
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
%   hazardORcentroids:  a TR hazard event set, or a filename of a saved
%                       one. The centroids for the ls hazard set will be
%                       constructed from the lat lons.
%                    OR a centroids struct, which will also form the basis
%                       for the rainfall hazard set.
%                       > prompted for if not given
%                       NOTE that for stable performance, centroids should
%                       form a regular grid.
%   focus_areas:        struct array containing shape files defining focus
%                       areas of interest. Useful considering the
%                       limitation imposed by the necessity for a regular
%                       grid of centroids. Land slides will only be
%                       calculated starting at centroids inside the regions
%                       defined by these polygons. Can also be the pathname
%                       of a .shp or .mat file.
%   hazard_set_file: the name of the hazard set file
%       > prompted for if not given
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
%   centroidsORhazard:  if rainfall hazard set is provided as input arg,
%                       then output is the centroids struct which is
%                       constructed from this hazard set, and to which
%                       topographic parameters have been assigned. If
%                       centroids are given as input, then output is the
%                       rainfall hazard set which has been calculated for
%                       these centroids.
% MODIFICATION HISTORY:
%   Gilles Stassen, gillesstassen@hotmail.com, 20150330
%   Gilles Stassen, 20150708, renamed MS->LS, time-dependent soil moisture, linear decay according to ET_mm_day, separate input args combined into hazardORcentroids
%   Gilles Stassen, 20150710, factor of safety calculation vectorised over events
%   Gilles Stassen, 20150713, focus_areas added as input
% -

hazard = []; centroidsORhazard = [];    % init

global climada_global
if ~climada_init_vars, return; end

% check arguments
if ~exist('hazardORcentroids',  'var'),     hazardORcentroids   = [];   end
if ~exist('focus_areas',        'var'),     focus_areas         = [];   end
if ~exist('hazard_set_file',    'var'),     hazard_set_file     = [];   end
if ~exist('check_plots',        'var'),     check_plots         = 1;    end
if ~exist('chrono_check',       'var'),     chrono_check        = 1;    end
if ~exist('bool_hazard_check',  'var'),     bool_hazard_check   = 0;    end

module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];

if isempty(hazardORcentroids) || ~isstruct(hazardORcentroids)
    res = input('would you like to auto-generate a rainfall hazard set (a); or browse files (b)?\t','s');

    if strcmp(res,'a') 
        [admin_name, admin_shapes, ~, admin_shape_ndx,country_name] = climada_admin_name;
        % centroids = climada_generate_centroids([],0.09,0); % hard-wired 90m resolution for srtm DEM
        centroids = climada_90m_DEM(admin_shapes(admin_shape_ndx),'DL','NO_SAVE');
        
        [country_name, country_ISO3] = climada_country_name(country_name);
        centroids.country_name = repmat({country_name},size(centroids.centroid_ID));
        centroids.admin0_ISO3 = country_ISO3;   centroids.admin0_name = country_name;
        hazard_rf = climada_rf_hazard_set(centroids,[],'NO_SAVE',1);
    elseif strcmp(res,'b') 
        % prompt for RF hazard event set if not given
        hazard_rf_file=[module_data_dir filesep 'hazards' filesep '*.mat'];
        [filename, pathname] = uigetfile(hazard_rf_file,...
            'Select a rainfall hazard event set:');
        if isequal(filename,0) || isequal(pathname,0)
            return; % cancel
        else
            hazard_rf_file =fullfile(pathname,filename);
        end
        load(hazard_rf_file);
        hazard_rf = hazard; clear hazard
        % construct centroids from RF hazard
        centroids.lon = hazard_rf.lon;  
        centroids.lat = hazard_rf.lat;
        centroids.centroid_ID = hazard_rf.centroid_ID;
        centroids.comment = ['centroids created from: ' hazard_rf.comment];
    end
elseif isfield(hazardORcentroids,'intensity') && isfield(hazardORcentroids,'peril_ID') % input is hazard set
    % input is 
    hazard_rf = hazardORcentroids; clear hazardORcentroids
    % construct centroids from RF hazard
    centroids.lon = hazard_rf.lon;  
    centroids.lat = hazard_rf.lat;
    centroids.centroid_ID = hazard_rf.centroid_ID;
    centroids.comment = ['centroids created from: ' hazard_rf.comment];
else % input is centroids (probably) so auto-generate rainfall hazard
    centroids = hazardORcentroids; clear hazardORcentroids
    hazard_rf = climada_rf_hazard_set(centroids,[],'NO_SAVE',1);
    centroidsORhazard = hazard_rf;
end

if ~isempty(focus_areas)
    if ischar(focus_areas) && exist(focus_areas,'file') % pathname given as input
        [fP, fN, fE] = fileparts(focus_areas);
        if strcmp(fE,'.shp')
            focus_areas = climada_shaperead(focus_areas,1,0);
        elseif strcmp(fE,'.mat')
            tmp = load(focus_areas);
            focus_areas = tmp.(getelements(fieldnames(tmp),1));
            clear tmp
        else
            cprintf([1 0.5 0],'WARNING: invalid filetype %s - computation continues for all centroids\n',fE)
            focus_areas = [];
        end
    elseif ~(isstruct(focus_areas) && isfield(focus_areas,'X') && isfield(focus_areas,'Y'))
        cprintf([1 0.5 0],'WARNING: invalid input - focus areas must be a .shp or .mat file, and must have fields .X and .Y\n')
        cprintf([1 0.5 0],'\t\tcomputation continues for all centroids\n')
        focus_areas = [];
    end
end     


% for testing purposes
% snapshot_year = [2004 2005]; % 2006 2007];
% hazard_rf = climada_hazard_extract_event(hazard_rf,hazard_rf.event_ID(ismember(hazard_rf.yyyy, snapshot_year)));
% hazard_rf = climada_hazard_extract_event(hazard_rf,-[1:10]);

hazard                  =   hazard_rf;
hazard.peril_ID         =   'LS';
hazard.date             =   datestr(now);
hazard.filename         =   hazard_set_file;
hazard.matrix_density   =   0.01; % estimate
hazard.comment          =   '';
% allocate the hazard array (sparse, to manage memory)
% hazard.intensity = spalloc(hazard_rf.event_count,length(hazard_rf.centroid_ID),...
%     ceil(hazard.event_count*length(hazard.lon)*hazard.matrix_density));
hazard.intensity        =   zeros(size(hazard_rf.intensity));
hazard.units            =   'm';

hazard.source           =   zeros(size(hazard_rf.intensity));
hazard.deposit          =   zeros(size(hazard_rf.intensity));
hazard.slide_ID         =   zeros(size(hazard_rf.intensity));

if isfield(hazard,'rainfield_comment')
    hazard = rmfield(hazard, 'rainfield_comment');
end

centroids_n_flds = length(fieldnames(centroids));

% Step 1: Compute centroid elevation fro
if ~isfield(centroids,'elevation_m')
    centroids = climada_90m_DEM(centroids,'DL',[],[],0);
    hazard.elevation_m = centroids.elevation_m;
end

% Step 2: Calculate topographic wetness index (TWI) 
if ~isfield(centroids,'TWI')
    centroids = centroids_TWI(centroids, 0);
end

% Step 3: Delineate basins
if ~isfield(centroids,'basin_ID')
    if isfield(climada_global,'climada_global_ori')
        ori = climada_global;
        climada_global = climada_global.climada_global_ori;
    end
    centroids = centroids_basin_ID(centroids, 15, 0);
    if exist('ori','var'), climada_global = ori; end
end

% Step 3: Compute daily evapotranspiration (ET)
if ~isfield(centroids,'ET_mm_day')
    centroids = centroids_ET(centroids, 0);
end

% % Step 4: Assign soil wetness index (SWI)
% if ~isfield(centroids,'SWI') || force_recalc
%     centroids = centroids_SWI(centroids, check_plots);
% end

% Step 5: Assign available water-holding capacity of the soil (WHC)
if ~isfield(centroids,'WHC_mm')
    centroids = centroids_WHC(centroids, 0);
end

% Step 6: Assign soil bulk density values (BD)
if ~isfield(centroids, 'BD_kg_m3')
    centroids = centroids_BD(centroids,0);
end

% Step 7: add soil depth
if ~isfield(centroids, 'SD_m')
    centroids = centroids_SD(centroids,0);
end

% Step 8: add leaf area index
if ~isfield(centroids, 'LAI')
    centroids = centroids_LAI(centroids,0);
end

% if ~isfield(centroids, 'EC_Pa')
    % centroids.EC_Pa  =   6400 * ones(size(centroids.centroid_ID)); %1900 * mean(centroids.LAI,1); %proxy
    centroids.EC_Pa  =   7250 * ones(size(centroids.centroid_ID));
% end

% auto save if updated (i.e. new fields added)
% if isfield(centroids,'filename') && length(fieldnames(centroids)) > centroids_n_flds
%     fprintf('autosaving centroids with additional fields to %s \n',centroids.filename)
%     save(centroids.filename,'centroids')
% end
clear centroids_n_flds


basin_IDs           =   unique(centroids.basin_ID);
n_basins            =   length(basin_IDs);
n_events            =   hazard_rf.event_count;
n_centroids         =   length(centroids.centroid_ID);

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

% init for plotting data
for event_i = 1: n_events
    slide_data(event_i).X = [];
    slide_data(event_i).Y = [];
    slide_data(event_i).Z = [];
    slide_data(event_i).C = [];
end

% for progress mgmt
format_str	= '%s';
t0          = clock;
mod_step    = 1;

% introduce time dependence into soil moisture calculation if hazard_rf
% allows it (i.e. events are in chronological order and real)
if all(unique(hazard_rf.datenum) == hazard_rf.datenum) && chrono_check
    chrono = 1;
else
    chrono = 0;
end

% soil moisture
soil_moisture = NaN(size(hazard_rf.intensity)); % init
ET_mm_day     = ones(size(centroids.centroid_ID)).*2;% ET_mm_day(isnan(ET_mm_day)) = min(ET_mm_day(~isnan(ET_mm_day)))/2; % proxy
for basin_i = 1: n_basins
    
    % index of centroids belonging to basin_i
    c_ndx       =   (centroids.basin_ID == basin_IDs(basin_i)) & (centroids.onLand==1);
    
    % fl_score_sum                =   sum(centroids.flood_score(c_ndx));
    wet_index_sum               =   sum(centroids.TWI(c_ndx & ~isnan(centroids.TWI)));
    
    for event_i = 1 : n_events
        
        if event_i>1 && chrono % since soil moisture has been initialised to zero for first event
            dt = hazard_rf.datenum(event_i) - hazard_rf.datenum(event_i-1);
            soil_moisture(event_i,c_ndx) = nanmax(soil_moisture(event_i-1,c_ndx) - ET_mm_day(c_ndx) .* dt,zeros(1,sum(c_ndx)));
        end
        
        % progress mgmt
        if mod(event_i,mod_step) ==0
            mod_step        = 100;
            t_elapsed       = etime(clock,t0)/((basin_i-1)*n_events+event_i);
            n_remaining     = (n_basins-basin_i+1)*n_events-event_i;
            t_projected_sec = t_elapsed*n_remaining;
            if t_projected_sec<60
                msgstr = sprintf('computing soil moisture, est. %3.0f sec left (%i/%i basins): (%i/%i events)',t_projected_sec, basin_i, n_basins,event_i,n_events);
            else
                msgstr = sprintf('computing soil moisture, est. %3.1f min left (%i/%i basins): (%i/%i events)',t_projected_sec/60, basin_i, n_basins,event_i,n_events);
            end
            fprintf(format_str,msgstr);
            format_str = [repmat('\b',1,length(msgstr)) '%s'];
        end
        
        % index of rained-on centroids
        r_ndx   = hazard_rf.intensity(event_i,:) > 0;
        
        if ~any(r_ndx & c_ndx)
            continue
        end
        
        % index of potentially moist centroids, i.e. centroids lower in
        % elevation than highest r_ndx centroid
        ls_ndx    =     c_ndx & (centroids.elevation_m <= max(centroids.elevation_m(r_ndx & c_ndx)));
        ls_ndx    =     ls_ndx & ~isnan(centroids.TWI);
        
        rain_sum  =   sum(hazard_rf.intensity(event_i,r_ndx & c_ndx),2);
        
        % calculate incoming rainfall per centroid according to TWI
        Q_in = zeros(1,n_centroids); %init
        if wet_index_sum ~=0
            Q_in(ls_ndx) = rain_sum .*(centroids.TWI(ls_ndx) ./ wet_index_sum);
        else
            Q_in(ls_ndx) = rain_sum / sum(ls_ndx);
        end
        
        if chrono
            % add soil moisture to sm_init, cap at water holding capacity
            soil_moisture(event_i,c_ndx) = min(soil_moisture(event_i,c_ndx)+Q_in(c_ndx),centroids.WHC_mm(c_ndx));
        else
            % subtract evapotranspiration (assume events of length 1 day), cap at water holding capacity
            soil_moisture(event_i,c_ndx) = max(min(soil_moisture(event_i,c_ndx)-ET_mm_day(c_ndx),centroids.WHC_mm(c_ndx)),zeros(1,sum(c_ndx)));
        end
    end
end

fprintf(format_str,...
    sprintf('processing soil moisture at %i centroids for %i rainfall events took %i minutes\n',...
    n_centroids, n_events,round(etime(clock,t0)/60)));

% sort by elevation (work down the mountain)
[~,elev_ndx] = sort(centroids.elevation_m,'descend');

if ~isempty(focus_areas)
    % filter for centroids in focus area
    in = zeros(size(centroids.centroid_ID)); %init
    for a_i = 1:length(focus_areas)
        in = in | inpolygon(centroids.lon,centroids.lat,focus_areas(a_i).X,focus_areas(a_i).Y);
    end
    elev_ndx = elev_ndx(ismember(elev_ndx,find(in))); % only loop over desired centroids
    hazard.focus_areas      = zeros(size(centroids.centroid_ID));
    hazard.focus_areas(in)  = 1; % for identification later
end

% calculate factor of safety for each event
fprintf('calculating factor of safety for all events... ')
try % vectorised first, since faster, but takes a lot of memory for many centroids
    hazard.factor_of_safety = climada_ls_cell_failure(centroids, soil_moisture);
catch % so loop over events if matlab unable to handle
    for event_i =1:n_events
        hazard.factor_of_safety(event_i,:) = climada_ls_cell_failure(centroids,soil_moisture(event_i,:));
    end
end
fprintf('done\n')

% land slides
hazard.intensity        =   zeros(size(hazard_rf.intensity)); %init

% for progress mgmt
t0          = clock;
mod_step    = 1;
format_str  = '%s';

for event_i = 1 : n_events

    tmp_factor_of_safety = hazard.factor_of_safety(event_i,:);
    if ~any(tmp_factor_of_safety(~isnan(tmp_factor_of_safety)) < 1),  continue;   end
    % tmp_factor_of_safety(tmp_factor_of_safety < 1) = 0.0;
    
    % progress mgmt
    if mod(event_i,mod_step) ==0
        mod_step        = 10;
        t_elapsed       = etime(clock,t0)/event_i;
        n_remaining     = n_events-event_i;
        t_projected_sec = t_elapsed*n_remaining;
        if t_projected_sec<60
            msgstr = sprintf('generating land slides, est. %3.0f sec left (%i/%i events)',t_projected_sec, event_i,n_events);
        else
            msgstr = sprintf('generating land slides, est. %3.1f min left (%i/%i events)',t_projected_sec/60, event_i,n_events);
        end
        fprintf(format_str,msgstr);
        format_str = [repmat('\b',1,length(msgstr)) '%s'];
    end
    
    for centroid_i = elev_ndx

        % if factor of safety is less than rand, failure occurs at
        % centroid_i
        if tmp_factor_of_safety(1,centroid_i) < rand
            
            % set specific temp FoS to 1, to avoid recheck for other
            % clusters
            tmp_factor_of_safety(1,centroid_i) = 1;
            
            % index for centroids at lower elevation
            r = climada_geo_distance(centroids.lon(centroid_i),centroids.lat(centroid_i),...
                centroids.lon,centroids.lat);
            [~, r_ndx]  = sort(r,'ascend');
            min_r       = r(r_ndx(2));
            
            % initialise
            cluster_ndx  = zeros(size(centroids.centroid_ID));
            cluster_size = 0; cluster_ndx(centroid_i) = 1; i =0;
            sink_ndx     = centroids.centroid_ID == centroids.centroid_ID(centroid_i);
            
            while sum(cluster_ndx) > cluster_size % loop until cluster stops growing
                
                cluster_size    = sum(cluster_ndx);
                
                sink_ndx        = centroids.centroid_ID == centroids.sink_ID(sink_ndx);
                % cluster_ndx     = cluster_ndx |(sink_ndx & (tmp_factor_of_safety(event_i,:)./cluster_size) < rand);
                cluster_ndx     = cluster_ndx |(sink_ndx & (tmp_factor_of_safety(1,:) < rand));
                
                % set specific temp FoS to 1, to avoid recheck for other
                % clusters
                tmp_factor_of_safety(1,cluster_ndx) = 1;
                
            end % loop until cluster stops growing
            
            % index of last cluster cell for start of deposit
            cluster_end_ndx = sink_ndx;
            
            % ************************************************************
            % insist that at least 2 neighbouring cells must fail for a
            % land slide to occur
%             if cluster_size < 2
%                 continue;
%             end
            % ************************************************************
            
            % deposition of slide
            slide_ndx       = false(size(centroids.centroid_ID));
            slope           = centroids.slope_deg(centroid_i); %init
            sink_ndx        = centroids.centroid_ID == centroids.centroid_ID(centroid_i);
            i = 0;
            while slope > 2 ...                                     % slope is steep enough
                    && ~isempty(centroids.sink_ID(sink_ndx))...     % there exists a sink to deposit into
                    && ~isnan(centroids.sink_ID(sink_ndx)) ...      % the sink is defined
                    && i < (2*sum(cluster_ndx)) ...                 % limit deposition length to twice the source length
                    && slide_ndx(centroids.sink_ID(sink_ndx)) ==0   % there must not already be a land slide at the next sink
                
                % loop counter
                i = i+1;
                sink_ndx     = (centroids.centroid_ID == centroids.sink_ID(sink_ndx));  % find next sink
                slope        = centroids.slope_deg(sink_ndx);                           
                slide_ndx    = slide_ndx | sink_ndx;    % add sink to slide

%                 if     ~isempty(centroids.sink_ID(sink_ndx))...     % there exists a sink to deposit into
%                     && ~isnan(centroids.sink_ID(sink_ndx)) ...      % the sink is defined
%                     && i < (2*sum(cluster_ndx)) ...                 % limit deposition length to twice the source length
%                     && slide_ndx(centroids.sink_ID(sink_ndx)) ==0   % there must not already be a land slide at the next sink
%                     break;
%                 end
                
                % p = scatter(centroids.lon(deposit_ndx),centroids.lat(deposit_ndx),'filled');
                % set(p,'markerFaceColor','g','markerEdgeColor','g')
            end
            
            % set temp FoS to 1 for entire slide, to avoid recalculation
            tmp_factor_of_safety(1,cluster_ndx) = 1;
            
            if any(isnan(centroids.area_m2(slide_ndx)))
                centroids.area_m2(isnan(centroids.area_m2) & slide_ndx) = nanmean(centroids.area_m2(slide_ndx));
            end
            
            % for testing purposes
            deposit_ndx = slide_ndx & ~cluster_ndx;
            
            hazard.source (event_i,cluster_ndx) = centroids.area_m2(cluster_ndx) .* centroids.SD_m(cluster_ndx);
            hazard.deposit(event_i,deposit_ndx) = sum(hazard.intensity(event_i,cluster_ndx)) / sum(deposit_ndx);
            hazard.slide_ID(event_i,slide_ndx)  = str2num(sprintf('%i%i%i',basin_i,event_i,centroid_i));
            
            hazard.intensity(event_i,cluster_ndx) = - centroids.SD_m(cluster_ndx);
            
            % for better overview of the code
            lon = centroids.lon;
            lat = centroids.lat;
            
            % distances from end of cluster
            dep_dist_m = climada_geo_distance(lon(cluster_end_ndx),lat(cluster_end_ndx),...
                lon(deposit_ndx),lat(deposit_ndx));
            
            % deposit intensity factor depends on distance from end of cluster
            int_factor = 2.* (1- dep_dist_m./max(dep_dist_m));
            int_factor = 1; % or not..
            
            hazard.intensity(event_i,deposit_ndx) = int_factor .* ...
            (sum(hazard.source(event_i,cluster_ndx))/sum(centroids.area_m2(deposit_ndx)));
            
            X  = centroids.lon(slide_ndx);
            Y  = centroids.lat(slide_ndx);
            Z  = centroids.elevation_m(slide_ndx);
            C  = hazard.intensity(event_i,slide_ndx);
            
            % attempt to smooth out the land slides - doesn't work so
            % well...
            smooth_factor = 1;
            if smooth_factor >1
                X_ = []; Y_ = []; Z_ = []; C_ = []; S = 0; Sq =[];
                % parameterise length along slide
                for s_i = 2:sum(slide_ndx)
                    S(end+1) = S(end) + sqrt((X(s_i)-X(s_i-1))^2 + (Y(s_i) - Y(s_i-1))^2 + (Z(s_i) - Z(s_i-1))^2);
                end
                for s_i = 1:length(S)-1
                    Sq(end+1:end+smooth_factor) = linspace(S(s_i),S(s_i+1),smooth_factor);
                end
                X_ = interp1(S,X,Sq,'spline');
                Y_ = interp1(S,Y,Sq,'spline');
                Z_ = interp1(S,Z,Sq,'spline');
                C_ = interp1(S,C,Sq,'spline');
                
                X = X_; Y = Y_; Z = Z_; C = C_;clear X_ Y_ Z_ C_
            end
        
            % vectors for easy plotting
            slide_data(event_i).X = [slide_data(event_i).X NaN X];
            slide_data(event_i).Y = [slide_data(event_i).Y NaN Y];
            slide_data(event_i).Z = [slide_data(event_i).Z NaN Z];
            slide_data(event_i).C = [slide_data(event_i).C NaN C];
        end
    end
end
hazard.intensity         = sparse(hazard.intensity);
hazard.matrix_density    = nnz(hazard.intensity)/numel(hazard.intensity);

fprintf(format_str,...
    sprintf('generating land slides at %i centroids for %i rainfall events took %i minutes\n',...
    n_centroids, n_events,round(etime(clock,t0)/60)));

if bool_hazard_check % only ones (land slide) and zeroes (no land slide)
    hazard.intensity(hazard.intensity ~=0) = 1;
    hazard.unit = '';
end

% filter out all-zero events
nz_events_ndx = find(sum(abs(hazard.intensity),2) ~= 0);
if length(nz_events_ndx) < hazard.event_count
    fprintf('filtering hazard set for %i non-zero events... ',length(nz_events_ndx))
    hazard      = climada_hazard_extract_event(hazard,nz_events_ndx);
    slide_data  = slide_data(nz_events_ndx);
    fprintf('done\n')
end
hazard.frequency = ones(size(hazard.event_ID))./hazard.event_count;

[~, centroids_name, ~] = fileparts(centroids.filename);
if ~isfield(centroids,'admin0_name'), centroids.admin0_name = ''; end
hazard.comment = sprintf('LS hazard set %s, centroids: %s, created on: %s',centroids.admin0_name, centroids_name,datestr(now,'HH:MM ddmmyyyy'));

if exist('res','var') && nargout == 2
    % neither centroids nor RF hazard were provided as input, so return both
    % in cell array
    centroidsORhazard = {centroids; hazard_rf};
elseif isfield(centroids,'comment') && ...
        strcmp(centroids.comment,['centroids created from: ' hazard_rf.comment])
    % centroids created from hazard input
    centroidsORhazard = centroids;
end

% prompt for ms_hazard_save_file if not given
if isempty(hazard_set_file) % local GUI
    hazard_set_file=[module_data_dir filesep 'hazards' filesep 'LS_hazard.mat'];
    [filename, pathname] = uiputfile(hazard_set_file, 'Save new land slide hazard event set as:');
    if isequal(filename,0) || isequal(pathname,0)
        hazard_set_file = 'NO_SAVE'; % cancel
        cprintf([0 0 1],'NOTE: hazard not saved\n')
    else
        hazard_set_file=fullfile(pathname,filename);
    end
end

if ~isempty(hazard_set_file) && ~strcmp(hazard_set_file,'NO_SAVE')
    fprintf('saving LS hazard set as %s\n',hazard_set_file);
    hazard.filename         =   hazard_set_file;
    try
        save(hazard_set_file,'hazard')
    end
    file_dir = dir(hazard_set_file);
    if isempty(file_dir) || file_dir.datenum < now -1 % one day
        fprintf('attempting to save using ''-v7.3'' switch... ')
        try
            save(hazard_set_file,'hazard','-v7.3')
        end
        if ~isempty(file_dir) && ~(file_dir.datenum < now -1)
            fprintf('success\n')
        else
            fprintf('failed\n')
        end
    end
end

% Plots
if check_plots 
    figure('Name', 'LS hazard set (largest event)','color', 'w')
    
    % get largest event
    event_sum       =   sum(abs(hazard.intensity) >0,2);
    [~,sort_ndx]    =   sort(event_sum,'descend');
    max_event       =   sort_ndx(1);
    hazard_intensity=   full(hazard.intensity(max_event,:));
    
    % also plot DEM for land slide hazard
    % construct regular grid
    [x, y]  = meshgrid(unique(hazard.lon),unique(hazard.lat));
    z       = griddata(hazard.lon,hazard.lat,centroids.elevation_m,x,y,'linear');
    xy(:,:,1) = x;  xy(:,:,2) = y; lonlat(:,1) = hazard.lon; lonlat(:,2) = hazard.lat;
    for y_i = 1:size(y,1); 
        nanpts(y_i,:) = ~ismember(squeeze(xy(y_i,:,:)),lonlat,'rows');
    end
    z(nanpts) = NaN;
    [C,h]   = contourf(unique(hazard.lon),unique(hazard.lat),z,10);
    l_h     = clabel(C,h);
    for i=1:length(l_h)
        s = get(l_h(i),'String');
        s = str2num(s);
        s = sprintf('%4.1f',s);
        set(l_h(i),'String',s);
    end
    colormap(flipud(bone(50)))
    freezeColors
    hold on
    
    if isfield(centroids,'rivers')
        plot3(centroids.rivers.X,centroids.rivers.Y,ones(size(centroids.rivers.X)).*5,'b')
    end
    freezeColors
    
    % use surface plot to show varying colour along landslide
    x_ = slide_data(max_event).X;
    y_ = slide_data(max_event).Y;
    c_ = slide_data(max_event).C; % (slide_data(max_event).C + min(slide_data(max_event).C(~isnan(slide_data(max_event).C))))./max(slide_data(max_event).C(~isnan(slide_data(max_event).C)));
    z_ = ones(size(x_));
    s  = surface([x_;x_],[y_;y_],[z_;z_],[c_;c_],'edgecol','interp','linew',5, 'marker','o','markersize',1);
    colormap(climada_colormap('MS'))
    %colormap(hot)
    caxis([min(hazard_intensity) max(hazard_intensity)])
    cb = colorbar;
    ylabel(cb,sprintf('landslide depth [%s]',hazard.units))
    
   	title_str=sprintf('%s event %i on %s',hazard.peril_ID,max_event,datestr(hazard.datenum(max_event),'dddd dd mmmm yyyy'));
    title(title_str)
    xlabel('Longitude')
    ylabel('Latitude')
    hold off
    
    figure('Name', 'LS hazard set (largest event) 3D','color', 'w')
    surf(unique(hazard.lon),unique(hazard.lat),z./(111.12 * 1000));
    %colormap(landcolor)
    colormap(flipud(bone(50)))
    shading interp
    freezeColors
    hold on
    
    % use surface plot to show varying colour along landslide
    x_ = slide_data(max_event).X;
    y_ = slide_data(max_event).Y;
    z_ = slide_data(max_event).Z ./(111.12 * 1000); % convert to degrees to use equal axis
    c_ = slide_data(max_event).C;
    s  = surface([x_;x_],[y_;y_],[z_;z_],[c_;c_],'edgecol','interp','linew',5, 'marker','o','markersize',1);
    colormap(climada_colormap('MS'))
    %colormap(flipud(hot))
    caxis([min(hazard_intensity) max(hazard_intensity)])
    cb = colorbar;
    ylabel(cb,sprintf('landslide depth [%s]',hazard.units))
    
   	title_str=sprintf('%s event %i on %s',hazard.peril_ID,max_event,datestr(hazard.datenum(max_event),'dddd dd mmmm yyyy'));
    title(title_str)
    xlabel('Longitude')
    ylabel('Latitude')
    
    set(gca,'ztick',[], 'zcolor', 'w');
    view(60,30)
    axis equal
    grid off
    box off
    
end

return         
   
% Some Notes
% ==================
% *     Might be interesting to consider capping the number of landslides
%       in a certain location to (centroid resolution/typical land slide
%       width). This could also include some time dependence, i.e. land
%       slide can only occur in the same place again after X years.


% ==================
% the following is an old method for calculating land slides which may have
% a width greater than the centroids resolution
% ==================

% % loop counter
% i = i+1;
% 
% % ring of nearest, lower neighbours
% ring_ndx = (centroids.elevation_m < centroids.elevation_m(centroid_i) ...
%     & r < i*(min_r*1.4) & r > i*(min_r*0.8));
% 
% % aspect criterion
% aspect_diff = mod(centroids.aspect_deg - centroids.aspect_deg(centroid_i),360);
% aspect_ndx = cosd(aspect_diff) > 0;
% 
% % grow cluster to include failing cells in ring
% cluster_ndx = cluster_ndx |(aspect_ndx & ring_ndx & hazard.factor_of_safety(event_i,:) < rand);
% %cluster_ndx = logical(cluster_ndx);