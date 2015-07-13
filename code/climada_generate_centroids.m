function centroids = climada_generate_centroids(centroids_rectORcountry_nameORshapes, resolution_km, buffer_check, save_file, check_plot)
% climada generate high resolution centroids
% MODULE:
%   barisal_demo
% NAME:
%   climada_generate_HR_centroids
% PURPOSE:
%   Given a rectangle defining the location of interest, generate an evenly
%   spaced rectilinear grid of hazard centroids
% CALLING SEQUENCE:
%   centroids = climada_generate_centroids(centroids_rectORcountry_nameORshapes, resolution_km, buffer_check, save_file, check_plot)
% EXAMPLE:
%   centroids = climada_generate_centroids
%   centroids = climada_generate_centroids(centroids_rect, 0.5, 1, 'NO_SAVE', 1)
%   centroids = climada_generate_centroids('Netherlands', 1.0, 1, '', 0)
%   centroids = climada_generate_centroids(shapes, 1.0, 1, '', 0)
% INPUTS:
%   centroids_rectORcountry_nameORshapes [prompted for if not given] 
%   can be any one of the following:
%       centroids_rect: 4-element row vector defining the longitude and latitude
%                       limits of the study region
%       country_name:   the name of a country of interest, or ISO3 code
%       shapes:         any generic shapes struct, with fields .X and .Y
%                       defining lat and lon coords of study boundary 
%                       region respectively.
% OPTIONAL INPUT PARAMETERS:
%   resolution_km:  specify the centroid resolution (default = 1 km)
%   buffer_check:   specifies whether a lower resolution grid of centroids 
%                   is generated outside the boundary defined by shapes,
%                   or if the high resolution centroids fill the entire
%                   boudning box (default = 1). Can also be set to -1, for
%                   centroids that only fill the shape
%   save_file:      full pathname of save location. If set to 'AUTO',
%                   centroids are automatically saved in the climada global 
%                   data directory. If set to 'NO_SAVE', centroids will not 
%                   be saved. (default = 'AUTO')
% OUTPUTS:
%   centroids:      climada centroids struct with fields
%                     .Longitude
%                     .Latitude
%                     .onLand
%                     .centroid_ID
%                     .countryname - cell array same size as .centroid_ID
%                     .admin0_name - country name char array
%                     .admin0_ISO3 - ISO 3 country code
%                     .comment
% MODIFICATION HISTORY:
% Gilles Stassen, gillesstassen@hotmail.com, 20150119
% Gilles Stassen, 20150128, added .comment field
% Gilles Stassen, 20150326, added buffer
% Gilles Stassen, 20150408, increased buffer size
% Gilles Stassen, 20150416, save functionality, general shape input
% Gilles Stassen, 20150423, documentation updated
% Gilles Stassen, 20150703, struct array shape input
%-
centroids = [];

global climada_global
if ~climada_init_vars, return; end

if exist(climada_global.map_border_file,'file')
    load(climada_global.map_border_file)
    shapes_check = 1;
else
    shapes_check = 0;
end

country_check = 0;
country_name = '' ; ISO3 = ''; % init
if ~exist('centroids_rectORcountry_nameORshapes','var') || isempty(centroids_rectORcountry_nameORshapes)
    if ~shapes_check
        % need shapes for country input
        fprintf('ERROR: no world border info found \n')
        return
    end
    
    if length(shapes) == 1 && isfield(shapes,'BoundingBox')
        centroids_rect =[shapes.BoundingBox(:,1)' shapes.BoundingBox(:,2)'];
    else
        if isfield(climada_global,'climada_global_ori')
            tmp_climada_global = climada_global;
            climada_global = climada_global.climada_global_ori;
            load(climada_global.map_border_file)
        end
        [country_name,ISO3,shape_index] = climada_country_name('Single');
        if exist('tmp_climada_global','var'), climada_global = tmp_climada_global; clear tmp_climada_global;  end
        shapes         = shapes(shape_index);
        if isempty(country_name), return; end % error message already printed in climada_country_name
        bb             = [min([shapes(:).X]) min([shapes(:).Y])
                          max([shapes(:).X]) max([shapes(:).Y])];
        centroids_rect = [bb(:,1)' bb(:,2)']; clear bb
        country_check  = 1;
    end
    
elseif ischar(centroids_rectORcountry_nameORshapes) 
    if ~shapes_check
        % need shapes for country input
        fprintf('ERROR: no world border info found \n')
        return
    end
    %input is country name
    country_name = centroids_rectORcountry_nameORshapes; clear centroids_rectORcountry_name
    if isfield(climada_global,'climada_global_ori')
        tmp_climada_global = climada_global;
        climada_global = climada_global.climada_global_ori;
        load(climada_global.map_border_file)
    end
    [country_name,ISO3,shape_index] = climada_country_name(country_name);
    if exist('tmp_climada_global','var'), climada_global = tmp_climada_global; clear tmp_climada_global;  end
    if isempty(country_name),   return;             end % error message already printed in climada_country_name
    shapes         = shapes(shape_index);
    bb             = [min([shapes(:).X]) min([shapes(:).Y])
                      max([shapes(:).X]) max([shapes(:).Y])];
    centroids_rect = [bb(:,1)' bb(:,2)']; clear bb
    country_check  = 1;
elseif isnumeric(centroids_rectORcountry_nameORshapes) && length(centroids_rectORcountry_nameORshapes) == 4
    % input is centroids rect
    centroids_rect = centroids_rectORcountry_nameORshapes; clear centroids_rectORcountry_name
    buffer_check   = 0;
elseif isstruct(centroids_rectORcountry_nameORshapes)
    %input is shapes
    shapes = centroids_rectORcountry_nameORshapes; clear centroids_rectORcountry_name
    if isfield(shapes,'X') && isfield(shapes,'Y')
        shapes_check = 1;
        bb             =   [min([shapes(:).X]) min([shapes(:).Y])
                            max([shapes(:).X]) max([shapes(:).Y])];
        centroids_rect = [bb(:,1)' bb(:,2)']; clear bb
    else
        fprintf('ERROR: invalid input argument, shapes must have fields ''X'' and ''Y'' \n')
        return
    end
else
    fprintf('ERROR: invalid input argument \n')
    return
end

if ~exist('resolution_km','var') || isempty(resolution_km)
    fprintf('WARNING: no resolution specified, resorting to default 1 km \n')
    resolution_km = 1;
end
resolution_ang = resolution_km / (111.12);

if ~exist('buffer_check',   'var'),     buffer_check    = 1;        end
if ~exist('save_file','var')||isempty(save_file),save_file = 'AUTO';end
if ~exist('check_plot',     'var'),     check_plot      = 0;        end

min_lon = centroids_rect(1)-2*resolution_ang;
max_lon = centroids_rect(2)+2*resolution_ang;
min_lat = centroids_rect(3)-2*resolution_ang;
max_lat = centroids_rect(4)+2*resolution_ang;

lon = [min_lon:resolution_ang:max_lon];
lat = [min_lat:resolution_ang:max_lat];

[lon, lat] = meshgrid(lon,lat);

centroids.lon = reshape(lon,1,numel(lon));
centroids.lat = reshape(lat,1,numel(lat));

% n_lon = round((max_lon - min_lon)/(resolution_ang)) + 1;
% n_lat = round((max_lat - min_lat)/(resolution_ang)) + 1;
% 
% % construct regular grid
% for i = 0 : n_lon - 1
%     ndx = i * n_lat;
%     centroids.lat(1,ndx + 1 : ndx + n_lat)= (1:n_lat) .* resolution_ang + min_lat;
%     centroids.lon(1,ndx + 1 : ndx + n_lat)= (n_lon - i) .* resolution_ang + min_lon;
% end

if shapes_check
    if buffer_check ==1
        % take only high resolution centroids for country
        in   = zeros(size(centroids.lon));
        for i = 1:length(shapes)
            in   = in | inpolygon(centroids.lon,centroids.lat,shapes(i).X,shapes(i).Y);
        end
        centroids.lon = centroids.lon(in);
        centroids.lat = centroids.lat(in);
        
        % generate coarser grid for buffer centroids (resolution reduced by
        % factor of 4)
        buffer_resolution_ang = resolution_ang * 4;
        min_lon = centroids_rect(1)-buffer_resolution_ang;
        max_lon = centroids_rect(2)+buffer_resolution_ang;
        min_lat = centroids_rect(3)-buffer_resolution_ang;
        max_lat = centroids_rect(4)+buffer_resolution_ang;

        lon = [min_lon:buffer_resolution_ang:max_lon];
        lat = [min_lat:buffer_resolution_ang:max_lat];

        [lon, lat] = meshgrid(lon,lat);

        buffer.lon = reshape(lon,1,numel(lon));
        buffer.lat = reshape(lat,1,numel(lat));
        
        % take only buffer centroids outside country border
        out = ones(size(buffer.lon));
        for i = 1:length(shapes)
            out     = out & ~inpolygon(buffer.lon,buffer.lat,shapes(i).X,shapes(i).Y);
        end
        buffer.lon = buffer.lon(out);
        buffer.lat = buffer.lat(out);
        
        % concatenate
        buffer_logical      = [zeros(size(centroids.lon)) ones(size(buffer.lon))];
        centroids.lon       = [centroids.lon buffer.lon];
        centroids.lat       = [centroids.lat buffer.lat];
    elseif buffer_check == -1
        % take only high resolution centroids for country
        in   = zeros(size(centroids.lon));
        for i = 1:length(shapes)
            in   = in | inpolygon(centroids.lon,centroids.lat,shapes(i).X,shapes(i).Y);
        end
        centroids.lon = centroids.lon(in);
        centroids.lat = centroids.lat(in);
    end
    
    if ~isempty(country_name) 
        for i = 1 : length(centroids.lon)
            centroids.country_name{i} = country_name;
        end
        centroids.admin0_ISO3 = ISO3;
    else
        if shapes_check && isfield(shapes,'NAME')
            for i = 1 : length(centroids.lon)
                centroids.country_name{i} = shapes.NAME;
            end
            centroids.admin0_name = shapes.NAME;
        end
        if isfield(shapes,'ADM0_A3')
            centroids.admin0_ISO3 = shapes.ADM0_A3;
        end
        if isfield(shapes,'NAME_0')
            for i = 1 : length(centroids.lon)
                centroids.country_name{i} = shapes(1).NAME_0;
            end
            centroids.admin0_name = shapes(1).NAME_0;
        end
        if isfield(shapes,'ISO')
            centroids.admin0_ISO3 = shapes(1).ISO;
        end
    end
end

centroids.centroid_ID = 1:length(centroids.lon);

% coastline is terrible.
centroids.onLand = ones(size(centroids.centroid_ID));
if exist('buffer_logical','var')
    centroids.onLand(buffer_logical ==1) = 0;
end
% 
% coastline = climada_coastline_read;
% 
% if ~isempty(coastline)
%     in = inpolygon(centroids.lon,centroids.lat,coastline.lon,coastline.lat);
%     centroids.onLand        = true(size(centroids.centroid_ID));
%     centroids.onLand(~in & buffer_logical)   = false;
% else
%     fprintf('WARNING: coastline info not found, .onLand set to NaN \n')
% end
centroids.comment = sprintf('%3.2f km resolution centroids, created on %s', resolution_km,datestr(now,'dd/mm/yyyy'));

if ischar(save_file) && ~strcmp(save_file,'NO_SAVE')
    if strcmp(save_file,'AUTO')
        if country_check
            save_file = [climada_global.data_dir filesep 'system' filesep 'centroids_' ISO3 '_' datestr(now,'ddmmyy') '.mat'];
        else
            save_file = [climada_global.data_dir filesep 'system' filesep 'centroids_' datestr(now,'ddmmyy') '.mat'];
        end
        fprintf('autosaving centroids as %s \n', save_file)
        centroids.filename = save_file;
        save(save_file, 'centroids')
    else
        fprintf('saving centroids as %s \n', save_file)
        centroids.filename = save_file;
        save(save_file, 'centroids')        
    end
end

if check_plot
    figure('name','Centroids','color','w')
    hold on
    climada_plot_world_borders
    scatter(centroids.lon(centroids.onLand ==1),centroids.lat(centroids.onLand == 1),'xr');
    scatter(centroids.lon(centroids.onLand ==0),centroids.lat(centroids.onLand == 0),'ob');
    if exist('country_name','var') && ~isempty(country_name)
        title(sprintf('Centroids for %s',country_name));
    end
    axis([min_lon max_lon min_lat max_lat])
    xlabel('Longitude')
    ylabel('Latitude')
end

return
