function centroids = climada_centroids_getAlti3d(rectangle,sourcepath,savepath)

% Generate a landslide hazard set.
% MODULE:
%   flood
% NAME:
%   climada_centroids_getAlti3d
% PURPOSE:
%   
% CALLING SEQUENCE:
%   centroids = climada_centroids_getAlti3d(rectangle)
% EXAMPLE:
%   
% INPUTS:
%   
% OPTIONAL INPUT PARAMETERS:
%   
% OUTPUTS:
%   
% MODIFICATION HISTORY:
% Thomas Rölli, thomasroelli@gmail.com, 20180409, init 

% check arguments
if ~exist('rectangle', 'var'), rectangle = []; end
if ~exist('sourcepath', 'var'), sourcepath = []; end

rectangle = [8.2456-.05 8.2456+.05 46.8961-.05 46.8961+.05];

deg_km = 111.12; %length of 1 degree on Earth

srtm1corner = [8.1958333333333329,46.846111111111114]; 
srtm1Res = 1/3600; 

path = 'C:\Users\Simon Rölli\Desktop\climada\climada_data\alti3d\alti3d_WGS84.tif';
alti3D = GRIDobj(path);
[lon,lat] = getcoordinates(alti3D);

[lon,lat] = meshgrid(lon,lat);
lat = flipud(lat);
alti3D.Z = flipud(alti3D.Z);

dlat = abs(min(diff(lat(:,1)))); 
dlon = abs(min(diff(lon(1,:))));

dy = dlat*(deg_km * 1000);
dx = dlon*cosd(mean(lat(:,1)))*(deg_km * 1000); 

minlat2 = ceil(min(lat(:))/srtm1Res)*srtm1Res;
maxlat2 = floor(max(lat(:))/srtm1Res)*srtm1Res;
minlon2 = ceil(min(lon(:))/srtm1Res)*srtm1Res;
maxlon2 = floor(max(lon(:))/srtm1Res)*srtm1Res;

lat2 = minlat2:srtm1Res:maxlat2;
lon2 = minlon2:srtm1Res:maxlon2;

minlat3 = lat2(find(lat2>rectangle(3),1,'first'));
maxlat3 = lat2(find(lat2<rectangle(4),1,'last'));
minlon3 = lon2(find(lon2>rectangle(1),1,'first'));
maxlon3 = lon2(find(lon2<rectangle(2),1,'last'));

lat3 = (minlat3:srtm1Res:maxlat3)';
lon3 = minlon3:srtm1Res:maxlon3;

[lon3,lat3] = meshgrid(lon3,lat3);

dem_grid = interp2(lon,lat,alti3D.Z,lon3,lat3,'linear');


disp('end')
end