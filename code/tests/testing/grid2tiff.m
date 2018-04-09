function grid2tiff(lon,lat,grid,filename)
%not a climada function
%transforms gridded data(n_lon,n_lat) into geotiff, which than can be
%used in ArcGIS
%intput:
%    lon,lat,grid:  grid(n_lon,n_lat) which should be transformed in a
%                   tiff-file
%    filename:      string of the name of the tiff file (without '.tif'). File is saved
%                   in '.\climada_data\geotiff'. If not given --> user is
%                   asked for path and filename.

global climada_global
if ~climada_init_vars, return; end

% check arguments
if ~exist('lon', 'var'), lon = []; end
if ~exist('lat', 'var'), lat = []; end
if ~exist('grid', 'var'), grid = []; end
if ~exist('filename', 'var'), filename = []; end

%check if folder to save geotiff exists
folder = [climada_global.data_dir filesep 'geotiff'];
if ~exist(folder)
    mkdir(folder);
end

%hazard_set_file      = [climada_global.data_dir filesep 'hazards' filesep 'LSXX_hazard.tif'];
%get path to save
if isempty(filename) % local GUI
    path      = [folder filesep 'grid2tiff.tif'];
    [filename,folder] = uiputfile(path, 'Save geotiff as:');
    if isequal(filename,0) || isequal(folder,0)
        return; % cancel
    else
        savepth = fullfile(folder,filename);
    end
else
    savepth = [folder filesep filename '.tif'];
end


DEM = GRIDobj(lon,lat,grid);
GRIDobj2geotiff(DEM,savepth);

fprintf('%s.tif saved in: \n %s',filename,folder)



end
