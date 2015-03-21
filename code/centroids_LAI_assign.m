function centroids=centroids_LAI_assign(centroids, LAI_img_filename, check_plots, force_recalc)
% Assign Leaf Area Indices (LAIs) to given centroids
% MODULE:
%   flood
% NAME:
%	centroids_LAI_assign
% PURPOSE:
%   Determine Leaf Area Indices (LAIs) for given centroids by
%   reading a grayscale picture of global annual mean LAIs. 
%   Leaf area index (LAI) is a dimensionless quantity that characterizes 
%   plant canopies. It is defined as the one-sided green leaf area per unit 
%   ground surface area (LAI = leaf area / ground area, m^2/m^2).
%   For more information on the LAI see 
%       http://en.wikipedia.org/wiki/Leaf_area_index
%   The Moderate Resolution Imaging Spectroradiometer (MODIS) aboard NASA's 
%   Terra and Aqua satellites collects global LAI data on a daily basis. 
%   Values range from 0 to 7 square meters of leaf area per square meter of 
%   land surface.
%   Global monthly mean LAI images in grayscale can be downloaded here:
%       http://neo.sci.gsfc.nasa.gov/view.php?datasetId=MOD15A2_M_LAI
% NOTE: Since only monthly mean data are provided under the link above, a
% grayscale picture of the annual mean LAI needs to be produced manually by
% adding monthly pictures and taking the average. 
% 
% CALLING SEQUENCE:
%   centroids = centroids_LAI_assign(centroids, LAI_img_filename)
% EXAMPLE:
%   centroids = centroids_LAI_assign(centroids)
% INPUTS:
%  centroids: Climada centroids struct; the following fields are required:
%         .lat:           Latitude
%         .lon:           Longitude
%         .centroid_ID:   centroid ID
%         .admin0_name    Country name
% OPTIONAL INPUT PARAMETERS:
%   LAI_img_filename: File name of the global mean LAI image. If empty, a
%   default image is used, and if it does not exists, the function prompts
%   for the image filename. 
%   check_plots: whether a plot should be drawn (=1) or not (=0; default)
%   force_recalc: if set to 1, LAIs are calculated even if they
%   already exist (default is 0)
% OUTPUTS: 
%   centroids: The input centroids structure with an additional field
%   'LAI', which contains the Leaf Area Index for each centroid
% NOTE: 
% MODIFICATION HISTORY:
% Melanie Bieli, melanie.bieli@bluewin.ch, 20150318, initial


% import/setup global variables
global climada_global
if ~climada_init_vars,return;end;

% check arguments
if ~exist('centroids','var') || isempty(centroids)
    cprintf([206 50 50]/255,['Error: Missing input ''centroids'' to '...
        'function centroids_LAI_assign, can''t proceed.\n']);
    return;
end
if ~exist('LAI_img_filename','var'),LAI_img_filename=''; end
if ~exist('check_plots', 'var') || isempty(check_plots) 
    check_plots = 0;    end
if ~exist('force_recalc','var') || isempty(force_recalc), 
    force_recalc = 0;  end

% locate the module's data folder
module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];

% PARAMETERS
%
% the file with the LAI data
default_LAI_img_filename=[module_data_dir filesep 'system' ...
    filesep 'LAI_global_average_2014.tif'];
min_South=-90; % degree, defined on the webpage above
max_North= 90; % defined on the webpage above
%
% prepare bounding box (for speedup, we only want to look at a section of
% the global LAI map
% bbox=[minlon minlat maxlon maxlat]
bbox = [min(centroids.lon), min(centroids.lat),...
    max(centroids.lon), max(centroids.lat)];
%-


% We only calculate the LAIs if the centroids do not come equipped
% with them or if force_recalc is set to 1
if ~isfield(centroids,'LAI') || force_recalc
    
    % Find LAI picture if not specified as an input argument
    % Try a default image filename first and then prompt for the file if the
    % default file does not exist
    if isempty(LAI_img_filename)
        [fP,fN,~] = fileparts(default_LAI_img_filename);
        default_LAI_img_filename_mat = [fP filesep fN '.mat'];
        if climada_check_matfile(default_LAI_img_filename_mat,'file')
            load(default_LAI_img_filename_mat);
        elseif exist(default_LAI_img_filename,'file')
            LAI_img_filename = default_LAI_img_filename;
        else
            % Prompt for image file
            LAI_img_filename=[module_data_dir filesep '*.tif'];
            [filename, pathname] = uigetfile(LAI_img_filename, 'Select LAI image:');
            if isequal(filename,0) || isequal(pathname,0) % Cancel pressed
                return
            else
                LAI_img_filename=fullfile(pathname,filename);
            end
        end
    end
    
    % Read the LAI image
    if ~isempty(LAI_img_filename)
        full_img=imread(LAI_img_filename);
        full_img=full_img(end:-1:1,:); % switch for correct order in latitude (images are saved 'upside down')
        % save the LAI image if no .mat file of it exists yet
        [fP,fN]=fileparts(LAI_img_filename);
        LAI_img_filename_mat=[fP filesep fN '.mat'];
        save(LAI_img_filename_mat,'full_img'); % for fast access next time
    end     
    
    % White means no data available, hence we set all white pixels to 0
    % (i.e., to black)
    full_img(full_img==255) = 0;
    
    % Crop the image to the centroids' bounding box (for speedup)
    xx = 360*(1:size(full_img,2))/size(full_img,2)+(-180); % -180..180
    yy = (max_North-min_South)*(1:size(full_img,1))/size(full_img,1)+min_South;
    pos_x = find(xx>=bbox(1) & xx<=bbox(3));
    pos_y = find(yy>=bbox(2) & yy<=bbox(4));
    img = full_img(pos_y,pos_x);
    xx = xx(pos_x);
    yy = yy(pos_y);
    [X,Y] = meshgrid(xx,yy); % construct regular grid
    
    % determine the pixel closest to each centroid
    n_centroids = length(centroids.centroid_ID);
    X_1D = reshape(X,[1,numel(X)]);
    Y_1D = reshape(Y,[1,numel(Y)]);
    fprintf('Assigning Leaf Area Indices to centroids...\n')
    for centroid_i=1:n_centroids
        distances=climada_geo_distance(centroids.lon(centroid_i),...
            centroids.lat(centroid_i),X_1D,Y_1D);
        [~,min_dist_index] = min(distances);
        [img_row, img_col] = ind2sub(size(X),min_dist_index);
        centroids.LAI(centroid_i) = img(img_row, img_col);
    end
    
    % convert range of greyscale values (0 to 255) into actual Leaf Area
    % Indices (LAIs) ranging from 0 to 7
    centroids.LAI = centroids.LAI.*(7/255);
    
    if check_plots
        % convert to double (from uint8)
        LAIs=double(img);
        % plot the image (kind of 'georeferenced')
        pcolor(X,Y,LAIs.*(7/255)); colorbar
        shading flat
        if isfield(centroids,'admin0_name')
            title_string = sprintf('Leaf Area Index (LAI), %s', ...
                centroids.admin0_name);
        else title_string = 'Leaf Area Index (LAI)';
        end
        title(title_string)
        hold on
        set(gcf,'Color',[1 1 1]) % white figure background
    end % if check_plots
    
else
    % centroids already have a field 'LAI'
    cprintf([23 158 58]/255,'Skipped - centroids already have LAIs.\n')
end % if isfield(centroids,'LAI')
    
end % centroids_LAI_assign




