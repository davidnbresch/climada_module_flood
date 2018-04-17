function centroids = climada_centroids_getAlti3d(rectangle,resolution,sourcepath,interpol)

% Get elevation data from alti3d for desired region
% MODULE:
%   flood
% NAME:
%   climada_centroids_getAlti3d
% PURPOSE:
%   It derives elevation data from Alti3d by interpolation. The new lon/lat
%   grid is such that it fits the centroid-coordinates of the SRTM1 or
%   SRTM3 data.
% CALLING SEQUENCE:
%   centroids = climada_centroids_getAlti3d(rectangle,sourcepath,savepath)
% EXAMPLE:
%   resolution of 1 arcsec --> SRTM1 ca. 30m (latitudinal)
%   centroids = climada_centroids_getAlti3d([8.2456-.05 8.2456+.05
%       46.8961-.05 46.8961+.05],1/3600)
%   Resolution of 1/3arcsec --> ca. 10m (latitudinal)
%   centroids = climada_centroids_getAlti3d([8.2456-.05 8.2456+.05 46.8961-.05 46.8961+.05],1/(3*3600))
%   larger study area
%   centroids = climada_centroids_getAlti3d([8.111086 8.341411 46.815369 46.946240],1/(3*3600))
% INPUTS:
%   rectangle:      rectangle[lonmin,lonmax,latmin,latmax] which defines
%                   regions for which centroid-grid should be created within.
%   resolution:     lat/lon resolution in arcsec. Should be a factor of
%                   1 arcsec to be compareable to other resolutions (eg.
%                   SRTM1 or SRTM3)
% OPTIONAL INPUT PARAMETERS:
%   sourcepath:     Path of alti3d.tiff file. If not given user is ask to
%                   provide by a UI.
%   interpol:       String which defines interpolation method when deriving
%                   new elevation values of new coordinates from alti3d.
%                   Possible input: 'linear' (default), 'nearest', 'cubic',
%                   'makima', 'spline'.
% OUTPUTS:
%   centroids:      .lat --> coordinates in latitudinal direction such that
%                   modulo of coordinate and resolution is zero
%                   .lon --> coordinates in longitudinal direction such that
%                   modulo of coordinate and resolution is zero
%                   .centroids_ID
%                   .comment:'interpolated from ALTI3D
%                   .sourcefile --> path where alti3d.tif is saved
% MODIFICATION HISTORY:
% Thomas Rölli, thomasroelli@gmail.com, 20180409, init
% Thomas Rölli, thomasroelli@gmail.com, 20180411, define sourcepath,
%   interpolation method

global climada_global
if ~climada_init_vars, return; end

% check arguments
if ~exist('rectangle', 'var'), rectangle = []; end
if ~exist('resolution', 'var'), resolution = []; end
if ~exist('sourcepath', 'var'), sourcepath = []; end
if ~exist('interpol', 'var'), interpol = []; end

if isempty(interpol), interpol = 'linear'; end

%init output
centroids = [];

%ask for sourcefile if not given
if isempty(sourcepath) % local GUI
    path      = [climada_global.data_dir];
    [filename,folder] = uigetfile('*.tif','load alti3d from:',path);
    if isequal(filename,0) || isequal(folder,0)
        return; % cancel
    else
        sourcepath = fullfile(folder,filename);
    end
else
    sourcepath = sourcepath;
end

deg_km = 111.12; %length of 1 degree on Earth

%read in tiff file
alti3D = GRIDobj(sourcepath);

%get coordinates of original tiff-file
[lon,lat] = getcoordinates(alti3D);
[lon,lat] = meshgrid(lon,lat);
lat = flipud(lat);
alti3D.Z = flipud(alti3D.Z);

%get resolution in lat and lon direction in meters, assume constant lat (no curvature, longituional resolution constant) 
dlat = abs(min(diff(lat(:,1)))); 
dlon = abs(min(diff(lon(1,:))));
dy = dlat*(deg_km * 1000);
dx = dlon*cosd(mean(lat(:,1)))*(deg_km * 1000); 

%get coordinate system which would be gotten when coordinates are a factor of
%resolution
minlat2 = ceil(min(lat(:))/resolution)*resolution; %now modulo of coordinates system and resolution = 0
maxlat2 = floor(max(lat(:))/resolution)*resolution;
minlon2 = ceil(min(lon(:))/resolution)*resolution;
maxlon2 = floor(max(lon(:))/resolution)*resolution;
lat2 = minlat2:resolution:maxlat2;
lon2 = minlon2:resolution:maxlon2;

%clip coordinats within rectangle
minlat3 = lat2(find(lat2>rectangle(3),1,'first'));
maxlat3 = lat2(find(lat2<rectangle(4),1,'last'));
minlon3 = lon2(find(lon2>rectangle(1),1,'first'));
maxlon3 = lon2(find(lon2<rectangle(2),1,'last'));
lat3 = (minlat3:resolution:maxlat3)';
lon3 = minlon3:resolution:maxlon3;
[lon3,lat3] = meshgrid(lon3,lat3);

%interpolate new elevation data with new coordinates from alti3d
dem_grid = interp2(lon,lat,alti3D.Z,lon3,lat3,interpol);

%resolution in meters of new elevation data
dlat3 = abs(min(diff(lat3(:,1)))); 
dlon3 = abs(min(diff(lon3(1,:))));
dy3 = dlat3*(deg_km * 1000);
dx3 = dlon3*cosd(mean(lat3(:,1)))*(deg_km * 1000); 

%write in centroids structure
centroids.lon = lon3(:)';
centroids.lat = lat3(:)';
centroids.elevation_m = dem_grid(:)';
centroids.centroid_ID = 1:numel(centroids.lon);
centroids.comment = 'Interpolated from ALTI3D';
centroids.sourcefile = sourcepath;

%disp('end')
end