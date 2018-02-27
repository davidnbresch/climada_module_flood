function mult_flow = climada_ls_multipleflow(centroids,hazard,exponent,test)

% Calculating of outflow proportion from each cell to its neighbours.
% Calculation according to Horton et al. (2013)
% MODULE:
%   flood
% NAME:
%   climada_ls_hazard_multipleflow
% PURPOSE:
%   
% CALLING SEQUENCE:
%   climada_ls_hazard_sets
% EXAMPLE:
%   
% INPUTS: 
%   centroids: a climada centroids stucture (including topographical
%              information)
% OPTIONAL INPUT PARAMETERS:
%   exponent:   variable exponent; is controlling the divergence of the flow
%               x=1: the spreading is similar to the multiple flow direction
%               x towards infinity: spreading similar to the single flow direction
%               Claessens et al. (2005) suggest x=4 for debris flow 
%   test:       If set true: a centroids structure with lon, lat and
%               elevation_m is constructed according to a test DEM (see
%               climada_ls_testDEM)
% OUTPUTS:
%   mult_flow:  8-D matrix with outflow proportion in each direction (each
%               neighbour-cell)
%  
% MODIFICATION HISTORY:
% Thomas Rölli, thomasroelli@gmail.com, 20180201, init
% Thomas Rölli, thomasroelli@gmail.com, 20180202, calculation of outflow
%   proportion
% Thomas Rölli, thomasroelli@gmail.com, 20180208, implementation of flow
%   path
% Thomas Rölli, thomasroelli@gmail.com, 20180214, implementaiton of
%   friction, outflow distance

%remove afterwards; load centroids and hazard
%load('C:\Users\Simon Rölli\Desktop\climada\climada_data\hazards\_LS_Sarnen_hazard.mat')
load('C:\Users\Simon Rölli\Desktop\climada\climada_data\hazards\_LS_Sarnen_srtm1_hazard.mat')

global climada_global
if ~climada_init_vars, return; end
% check arguments
if ~exist('centroids', 'var'), centroids = []; end
if ~exist('hazard', 'var'), hazard = []; end
if ~exist('exponent', 'var'), exponent = []; end
if ~exist('test', 'var'), test = []; end

% PARAMETERS 
if isempty(exponent); exponent = 4; end
if isempty(test); test = false; end

if test
   [centroids,hazard] = climada_ls_testDEM();
else
   %load('C:\Users\Simon Rölli\Desktop\climada\climada_data\centroids\_LS_Sarnen_centroids.mat')
   load('C:\Users\Simon Rölli\Desktop\climada\climada_data\centroids\_LS_Sarnen_srtm1_centroids.mat') 
end

%get dimension of grid field from lon/lat coordinates
%and reshap needed vectors --> easier to handel in grid format than in
%vector; only possible for regular placed gridpoints
%also flip matrix up down to have high latitude in first row of matrix
n_lon = numel(unique(centroids.lon));
n_lat = numel(unique(centroids.lat));
lon = flipud(reshape(centroids.lon,n_lat,n_lon));
lat = flipud(reshape(centroids.lat,n_lat,n_lon));
elevation = flipud(reshape(centroids.elevation_m,n_lat,n_lon));
intensity = logical(zeros(n_lat,n_lon,hazard.event_count));
for i = 1:hazard.event_count
    intensity(:,:,i) = flipud(reshape(hazard.intensity(i,:),n_lat,n_lon));
end

%calculate gradients from each cell to its 8 neighbours
[gradients,horDist,verDist] = climada_centroids_gradients(lon,lat,elevation);

% mutiplying by -1 --> outflow should be positive; neighbours with inflow are set to zero
gradients = gradients*-1;
gradients(gradients < 0) = 0; 


%%% calculate sum of all outflow cells
gradients_sum = sum(gradients.^exponent,3);
gradients_sum(gradients_sum==0) = 1; %prevent division by 0

%%% calculate multidirectional outflow proportion
% (tan(beta_i)^x/sum(tan(beta_i)^x(i= 1 to 8))
mult_flow = (gradients.^exponent)./gradients_sum;


%landslide path is calculated and the intensity spreaded according to the 
%multiple flow path and a simplified friction model
spread = climada_ls_spread(intensity,mult_flow,horDist,verDist,false);
spreadFri = climada_ls_spread(intensity,mult_flow,horDist,verDist,true);

%make some plots just to test some things
figure
surf(lon,lat,elevation,spread(:,:,1));
figure
surf(lon,lat,elevation,spreadFri(:,:,1));


%affected area: plot to see where intensity values are greater than 0
aff_area = zeros(n_lon,n_lat);
aff_area(spreadFri(:,:,1)>0) =1;
figure
surf(lon,lat,elevation,aff_area)
%surf(lon,lat,elevation,aff_area,'LineStyle','none')

%affected area: plot without friction
aff_area(spread(:,:,1)>0) =1;
figure
surf(lon,lat,elevation,aff_area)

disp('hier');













end
