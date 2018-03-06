function climada_ls_test_flowParameters()

% Read in generated landslide sets (with intensity coordinates (lat/lon) 
% and elevation. With this information try out flowpath with different 
% parameters. Different analyses for master thesis. Not for the enduser.
% MODULE:
%   
% NAME:
%   
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
% Thomas Rölli, thomasroelli@gmail.com, 20180306, init 

%load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards\_LS_Sarnen_hazard.mat')
load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_srtm1_hazard.mat')

%load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards\_LS_Sarnen_centroids.mat')
load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_srtm1_centroids.mat')

%get gridded datasets
n_lon = numel(unique(centroids.lon));
n_lat = numel(unique(centroids.lat));
lon = reshape(centroids.lon,n_lat,n_lon);
lat = reshape(centroids.lat,n_lat,n_lon);
elevation = reshape(centroids.elevation_m,n_lat,n_lon);
intensity = logical(zeros(n_lat,n_lon,hazard.event_count));
for i = 1:hazard.event_count
    intensity(:,:,i) = reshape(hazard.intensity(i,:),n_lat,n_lon);
end

%default parameters (exponent multipleflow = 25; max slide velocity = 8
%friction parameter phi = 18
% example: spreaded = climada_ls_flowpath(lon,lat,elevation,intensity,exponent,v_max,phi);
spreaded = climada_ls_flowpath(lon,lat,elevation,intensity);
figure
surface(lon,lat,elevation,spreaded(:,:,1))
figure
surface(lon,lat,elevation,double(spreaded(:,:,1)>0))


disp('hier')



end
