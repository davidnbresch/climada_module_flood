function hazard = climada_rf_hazard_set(precip_file, centroids, hazard_set_file, check_plots)
% climada RF hazard event set
% NAME:
%   climada_RF_hazard_set
% PURPOSE:
%   Construct hazard set structure from global historical daily precipitation data
%   ftp://ftp.cgd.ucar.edu/archive/PRECIP/GPCP_1DD_v1.2_199610-201407.nc.gz
%   http://www.esrl.noaa.gov/psd/data/gridded/data.gpcc.html
% CALLING SEQUENCE:
%   hazard=climada_rf_hazard_set(precip_file, centroids, hazard_set_file, check_plots)
% EXAMPLE:
%   hazard=climada_rf_hazard_set
%   hazard=climada_rf_hazard_set([],centroids,[],1)
% INPUTS:
%   hazard_set_file: the name of the newly created rainfall (RF) hazard
%       > promted for if not given
% OPTIONAL INPUT PARAMETERS:
%   check_plots: whether to show plots (default = 1)
% OUTPUTS:
%   hazard:     a hazard event set, see core climada doc
%               also written to a .mat file (see hazard_set_file)
% MODIFICATION HISTORY:
% Gilles Stassen, gillesstassen@hotmail.com, 20150316
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% poor man's version to check arguments
if ~exist('precip_file',        'var'), precip_file         ='DL';      end
if ~exist('hazard_set_file',    'var'), hazard_set_file     ='';        end
if ~exist('centroids',          'var'), centroids           =[];        end
if ~exist('check_plots',        'var'), check_plots         =0;         end
% PARAMETERS

module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];
if ~isdir(module_data_dir),mkdir(fileparts(module_data_dir),'data');end % create the data dir, should it not exist (no further checking)

if isempty(precip_file) || ~ischar(precip_file)
    precip_file  = 'GPCP_1DD_v1.2_199610-201407.nc';
    try
        precip_file = subdir([module_data_dir filesep precip_file]);
    catch
        precip_file = subdir([root_dir filesep precip_file]);
        precip_dir = fileparts(precip_file.name);
    end
    init_time = datenum(1990,01,01); % hard wired for GPCP
    
    if isempty(precip_file)
        warn_msg = sprintf(['WARNING: GPCP precipitation data file not found \n' ...
            '\t \t would you like to browse existing files (b) or download directly (d)?']);
        response = input(warn_msg,'s');
        switch response
            case 'd'       
                fprintf('downloading and unzipping global daily precipitation data from NASA''s GCPC data set... ')
                precip_URL = 'ftp://ftp.cgd.ucar.edu/archive/PRECIP/GPCP_1DD_v1.2_199610-201407.nc.gz';
                precip_file = gunzip(precip_URL,[module_data_dir filesep 'precip_data']);
                fprintf('done \n')
            case 'b'
                [fN, fP] = uigetfile('*.nc', 'Select precipitation data file');
                if isequal(fN,0) || isequal(fP,0),
                    cprintf([0.25 0.25 1],['NOTE: click '...
                        '<a href="ftp://ftp.cgd.ucar.edu/archive/PRECIP/GPCP_1DD_v1.2_199610-201407.nc.gz">'...
                        'here</a> to download global monthly precipitation data \n']);
                    return;
                end
                init_time   = input('specify start time yyyy/mm/dd of data set: ','s');
                init_time   = datenum(init_time,'yyyy/mm/dd');
                precip_file = [fP fN];
            otherwise
                cprintf([1 0 0], 'ERROR: invalid response \n')
                return
        end
    end
    
    fprintf('reading global monthly precipitation data... ')
    precip_nc_info = ncinfo(precip_file.name);
    for fld_i =  1 : numel(precip_nc_info.Variables)
        var_name = precip_nc_info.Variables(fld_i).Name;
        tmp_p_d.(var_name)=ncread(precip_file.name,var_name);
    end 
    fprintf('done \n')
else
    if exist(precip_file, 'file')
        [fP, fN, fE] = fileparts(precip_file);
        if strcmp(fE,'.mat')
            load(precip_file);
        elseif strcmp(fE, '.nc')
            fprintf('reading global monthly precipitation data... ')
            precip_nc_info = ncinfo(precip_file.name);
            for fld_i =  1 : numel(precip_nc_info.Variables)
                var_name = precip_nc_info.Variables(fld_i).Name;
                tmp_p_d.(var_name)=ncread(precip_file.name,var_name);
            end
            fprintf('done\n')
            init_time   = input('specify start time yyyy/mm/dd of data set: ','s');
            init_time   = datenum(init_time,'yyyy/mm/dd');
        else
            cprintf([1 0 0], 'ERROR: invalid file type - input file must be .mat or .nc \n')
            return
        end
    end
end

if ~isstruct(centroids)
    centroids = climada_centroids_load;
end

% prep the region we need (rectangular region encompassing the hazard centroids)
centroids_rect=[min(centroids.lon) max(centroids.lon)...
    min(centroids.lat) max(centroids.lat)];

% determine fields
flds = fieldnames(tmp_p_d);
for fld_i = 1: numel(flds)
    if ndims(tmp_p_d.(flds{fld_i})) == 3
        tmp_p_d.precip = tmp_p_d.(flds{fld_i});
    end
end

% restructure gridded data in array format
[p_data.precip, p_data.lon, p_data.lat]=climada_grid2array(permute(tmp_p_d.precip,[1 2 3]), tmp_p_d.lon, tmp_p_d.lat);
p_data.time = tmp_p_d.time + init_time;

% select relevant spatial region
s_ndx = (p_data.lon >= floor(centroids_rect(1)) & p_data.lon <= ceil(centroids_rect(2))) & ...
    (p_data.lat >= floor(centroids_rect(3)) & p_data.lat <= ceil(centroids_rect(4)));

p_data.lon      = p_data.lon(s_ndx);
p_data.lat      = p_data.lat(s_ndx);
p_data.precip   = p_data.precip(s_ndx,:);
clear s_ndx %tmp_p_d


% prompt for hazard_set_file if not given
if isempty(hazard_set_file) % local GUI
    hazard_set_file=[module_data_dir filesep 'hazards' filesep 'RF_hazard.mat'];
    [fN, fP] = uiputfile(hazard_set_file, 'Save new RF hazard event set as:');
    if isequal(fN,0) || isequal(fP,0)
        return; % cancel
    else
        hazard_set_file=fullfile(fP,fN);
    end
end

orig_years = year(p_data.time(end)) - year(p_data.time(1))+1;

% fill the hazard structure
hazard.reference_year   = climada_global.present_reference_year; % default for present hazard is normally 2010
hazard.lon              = centroids.lon;
hazard.lat              = centroids.lat;
hazard.centroid_ID      = centroids.centroid_ID;
if isfield(centroids,'elevation_m'),hazard.elevation_m=centroids.elevation_m;end
hazard.orig_years       = orig_years;

hazard.orig_event_count = length(p_data.time);
hazard.event_count      = length(p_data.time);
hazard.event_ID         = 1:hazard.event_count;
hazard.orig_event_flag  = ones(1,hazard.event_count);
hazard.yyyy             = str2num(datestr(p_data.time,'yyyy'))';
hazard.mm               = str2num(datestr(p_data.time,'mm'))';
hazard.dd               = str2num(datestr(p_data.time,'dd'))';
hazard.datenum          = p_data.time';

event_frequency         = 1/(orig_years);
hazard.frequency        = ones(1,hazard.event_count)*event_frequency;

hazard.peril_ID         = 'RF';
hazard.comment          =sprintf('RF hazard event set, generated %s',datestr(now));

%hazard.intensity        = spalloc(hazard.event_count,length(hazard.lon),ceil(hazard.event_count*length(hazard.lon)*0.3));
hazard.intensity        = zeros(hazard.event_count,numel(centroids.centroid_ID));
fprintf('processing RF precipitation at centroids for %i events...\n',hazard.event_count)
mod_step = 10; format_str = '%s'; t0 = clock;
[LON, LAT] = meshgrid(tmp_p_d.lon,tmp_p_d.lat);
% LON = unique(p_data.lon);
% LAT = unique(p_data.lat);
for event_i = 1:hazard.event_count
    if any(any(tmp_p_d.precip(:,:,event_i)))
        % hazard.intensity(event_i,arr_i)=qinterp2(LON, LAT,squeeze(double(precip_grid.precip(:,:,event_i))'),hazard.lon,hazard.lat);
        % hazard.intensity(event_i,arr_i)=interp2(LON, LAT, squeeze(double(precip_grid.precip(:,:,event_i))'),hazard.lon(arr_i),hazard.lat(arr_i));
        hazard.intensity(event_i,:)=interp2(LON, LAT, squeeze(tmp_p_d.precip(:,:,event_i)'),hazard.lon(:),hazard.lat(:));
    end
    % the progress management
    if mod(event_i,mod_step)==0
        mod_step          = 100;
        t_elapsed_event   = etime(clock,t0)/event_i;
        events_remaining  = hazard.event_count-event_i;
        t_projected_sec   = t_elapsed_event*events_remaining;
        if t_projected_sec<60
            msgstr = sprintf('est. %3.0f sec left (%i/%i events)',t_projected_sec,event_i,hazard.event_count);
        else
            msgstr = sprintf('est. %3.1f min left (%i/%i events)',t_projected_sec/60,event_i,hazard.event_count);
        end
        fprintf(format_str,msgstr); % write progress to stdout
        format_str=[repmat('\b',1,length(msgstr)) '%s']; % back to begin of line
    end
end
fprintf(strcat('\b',format_str),'done \n')
hazard.intensity    = sparse(hazard.intensity);
fprintf(format_str,sprintf('processing rainfall at %i centroids for %i events took %3.1f seconds \n', ...
    numel(hazard.centroid_ID),hazard.event_count,etime(clock,t0)));

if isfield(hazard,'filename'),hazard.filename_source=hazard.filename;end
hazard.filename=hazard_set_file;
hazard.date=datestr(now);
hazard.matrix_density=nnz(hazard.intensity)/numel(hazard.intensity);
hazard.units='mm'; % store the SI unit of the hazard intensity
if ~isfield(hazard,'orig_event_count') % fix a minor issue with some hazard sets
    if isfield(hazard,'orig_event_flag')
        fprintf('field hazard.orig_event_count inferred from hazard.orig_event_flag\n')
        hazard.orig_event_count=sum(hazard.orig_event_flag);
    else
        fprintf('WARNING: no field hazard.orig_event_flag\n')
    end
end

if ~strcmp(hazard_set_file,'NO_SAVE');
    fprintf('saving RF rainfall hazard set as %s\n',hazard_set_file);
    try
        save(hazard_set_file,'hazard');
    catch
        cprintf([1 0 0], 'ERROR: can not write to file, try saving RF hazard manually \n')
    end
end

if check_plots,figure('color','w'); climada_hazard_plot_hr(hazard,0);end % show max rainfall over ALL events

return
