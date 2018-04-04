function [slope,aspect,area,TWI] = climada_centroids_scores(lon,lat,dem,topoTWI)
% Calculate slope, aspect and topographic wetness index
% MODULE:
%   flood
% NAME:
%	climada_centroids_scores
% PURPOSE:
%   Calculate flood scores and topographic wetness indices for given
%   DEM, that are on a regular grid.
%   To calculate the flow accumulation a function from TopoToolbox is used
%   which implements the multiple flow according to Freeman T.G. (1991) --> flow to
%   several cells is possible, according to the gradients between central
%   cell ant its neigbouring cells. The dispersion can be controlled with
%   an exponent (set to 1.1 for Freeman T.G. (1991)).
%
%   The Topographic wetness index (TWI), which can be derived from the
%   calculated flow accumulation, then provides a cost-efficient
%   alternative to flood determination by conventional hydrodynamic models.
%   For more information on TWI see
%       Pourali, S.H. et al. (2014): Topography Wetness Index Application
%       in Flood-Risk-Based Land Use Planning;
%       doi:10.1007/s12061-014-9130-2
% CALLING SEQUENCE:
%   centroids = climada_centroids_TWI_calc(centroids, check_plots)
% EXAMPLE:
%   centroids = climada_centroids_TWI_calc(centroids,0)
% INPUTS:
%   lat:    Latitude in gridded format (nxm) --> coordinates of DEM raster
%   lon:    longitute in gridded format (nxm) --> coordinates of DEM raster
%   dem:    digital elevation model in gridded format (nxm)
%   topoTWI: 0 or 1(default). If set to zero topoToolbox is used to derive
%           the topographic wetness index (respectively, the upstream contributing
%           area
% OPTIONAL INPUT PARAMETERS:
%   check_plots: whether plots of topography (elevation), slope, aspect
%   angle, and flood scores should be drawn (=1) or not (=0; default)
%   force_recalc: if set to 1, flood scores are calculated even if they
%   already exist (default is 0)
% OUTPUTS:
%   centroids: centroids with additional fields:
%       .FL_score:   assigns a number for flow accumulation to each centroid, and
%       .TWI:        assigns a topographic wetness index to each centroid
%       .slope_deg:  slope of every centroid, in degree
%       .area_m2:    area of every centroid, in square meters
%       .aspect_deg: aspect of every centroid, in degree
%       .sink_ID:    sink of every centroid, links to centroid_ID
% MODIFICATION HISTORY:
% Thomas Rölli, thomasroelli@gmail.com, 20180404, init with
%  TopoToolbox-functions
global climada_global

% check input arguments
if ~climada_init_vars; return; end
if ~exist('lon', 'var'),   lon = []; end
if ~exist('lat', 'var'),   lat = []; end
if ~exist('dem', 'var'),   dem = []; end
if ~exist('topoTWI', 'var') topoTWI = []; end
if isempty(topoTWI), topoTWI = 1; end

% Weighting factor to be applied in the calculation of flow accumulation
% according to Freeman (1991)
weighting_factor = 1.1;

%init output
slope = zeros(size(lon));
aspect = zeros(size(lon));
area = zeros(size(lon));
TWI = zeros(size(lon));

%init DEMobj
DEM = GRIDobj(lon,flipud(lat),dem);

%%%%%%%%Normal Vector%%%%%%%

% Calculate mean difference between centroids in x and y direction
% We assume a constant distance between degrees of latitude (lat), i.e. we
% ignore the earth's slightly ellipsoid shape. The difference in
% longitude (lon) is calculated based on:
% length of 1 degree of lon = cosine(lat) * length of degree at equator
% lon_singleton = [min(centroids.lon):min(diff(unique(centroids.lon))):max(centroids.lon)];
% lat_singleton = [min(centroids.lat):min(diff(unique(centroids.lat))):max(centroids.lat)];
% 20180222 remove interpolation (because of slight change of coordinates

lon_singleton = lon(1,:);
lat_singleton = lat(:,1)';

% assume cos(lat) doesn't vary much within small study region and take cos
% of mean latitude
x = lon_singleton .* (cos(mean(mean(lat))*pi/180)* 111.12 * 1000); 
y = lat_singleton .* (111.12 * 1000);

factor_f = 1;

dx = diff(lon(1,1:2))*factor_f * (cos(mean(mean(lat))*pi/180)* 111.12 * 1000);
dy = diff(lat(1:2,1))*factor_f * (111.12 * 1000);
 
% Calculate gradients in x and y direction in order to derive normal
% vectors of the grid cells
vx = zeros([size(dem) 3]);
vy = zeros([size(dem) 3]);

[dzdx, dzdy] = gradient(dem, x, y);
vx(:,:,3)=dzdx.*dx;     vx(:,:,1)=dx;
vy(:,:,3)=dzdy.*dy;     vy(:,:,2)=dy;

normal_vector   =   cross(vx,vy);
area            =   sqrt(sum(normal_vector.^2,3));

%Normalize Vector
for ndx=1:3, normal_vector(:,:,ndx)=normal_vector(:,:,ndx)./area; end;


%%%%%%slope and aspect%%%%%%

% Calculate slope and aspect angle
slope               =   acosd(normal_vector(:,:,3));
aspect              =   atan2d(-1 .* normal_vector(:,:,2),normal_vector(:,:,1))+180;
aspect(slope ==0)   =   nan;

%%%%%%%%TWI%%%%%%%%%%

if topoTWI
    %%%%%%%%%get TWI using topoToolbox%%%%%
    %calculate flow accumulation by using 
    flowAcc = climada_ls_flowacc(DEM,1);
    %calculate upstream contributing area by multiplying by unit area of a
    %cell
    tmp_slope = slope + 0.1; % we don't want -inf values for wet_index
    TWI = log((1+flowAcc)*(dx*dy)./tand(tmp_slope));
else
    % Now we determine flow accumulation (which, in contrast to the local
    % parameters slope and aspect, can only be caluclated from the global
    % neighbourhood)
    % General concept: add cell inflows from higher adjacent cells, starting
    % from the the specified cell and working up to the watershed boundaries
    % The outflow of each cell is assigned to the neighbouring cell to which
    % the gradient is steepest

    % shift_matrix contains indices that will be used to access the
    % neighbouring cells
    % shift_matrix = [0 -1;-1 -1;-1 0;1 -1;0 1;1 1;1 0;-1 1];   % ori Melanie Bieli
    % shift_matrix = [0 -1;-1 -1;-1 0;-1 1;0 1;1 1;1 0;-1 1];   % Jacob Anz 280715
    shift_matrix = [0 -1;-1 -1;-1 0;-1 1;0 1;1 1;1 0;1 -1]; %Thomas Rölli 20180223 [c(8,:) and c(4,:) was the same before]
    [a, b]=size(dem);
    % gradients are stored in a 3D-matrix - gradients(:,:,1) refers to the
    % eastern neighbour, gradients(:,:,2) to the southeastern neighbour etc.
    gradients = zeros(a,b,8);

    % For southeastern, southwestern, northwestern and northeastern neighbours,
    % the gradient is calculated by taking the elevation difference of the two
    % cells and dividing it by sqrt(2)*grid_cell_size
    % For the time being, we use dx as a mean value for grid cell size; for
    % NOAA's global 1 arc-minute ETOPO1 data, the horizontal resolution at the
    % equator would be approximately 1.85 kilometers):
    % 070715 - diagonal from pythagoras
    for c = 2:2:8
        cell_diag = sqrt(dx^2+dy^2);
        gradients(:,:,c) = (circshift(dem,[shift_matrix(c,1) shift_matrix(c,2)])...
            -dem)/cell_diag;
    end
    % For eastern, western, southern, and northern neighbours, the gradient is
    % simply calculated by taking the elevation difference of the two cells and
    % dividing it by the cell size

    for c = 1:2:7
        if c == 1 || c == 5, distance=dx; else distance=dy; end
        gradients(:,:,c) = (circshift(dem,[shift_matrix(c,1) shift_matrix(c,2)])-dem)/distance;
    end

    % set gradients which are directed out of the study region to zero; at
    % boarder cells

    gradients(numel(lat(:,1)),:,[2 3 4]) = 0; %upper boarder
    gradients(1,:,[6 7 8]) = 0; %lower boarder
    gradients(:,1,[4 5 6]) = 0; %left boarder
    gradients(:,numel(lon(1,:)),[1 2 8]) = 0; %right boarder
    gradients = atand(gradients);%/pi*2;

    %----------------------
    % Prep for the calculation of flow accumulation (resulting in flood
    % scores), where the outflow of each cell will be distributed among its
    % neighbouring cells based on the steepness of the respective gradients
    % 
    % 1) Take the absolute value of negative gradients (i.e., the ones that
    % cause outflow) and raise them to the power of the selected weighting
    % factor.
    % Flow of positive gradients is set to zero.

    outflow_weighted = (gradients.*(-1*gradients<0)).^weighting_factor;

    % 2) Sum up gradients such that for each grid cell, outflow_gradients_sum
    % contains the sum of the gradients to the grid cell's 8 neighbouring cells
    outflow_gradients_sum = sum(outflow_weighted,3);
    outflow_gradients_sum(outflow_gradients_sum==0) = 1; %prevent division by 0

    % 3) Normalize the flow such that the neighboring cell gradients of each
    % direction (east, southeast, south, etc.) add up to 1. These fractions
    % will be used to distribute the flow to all neighbours of a cell
    for i = 1:8
        outflow_weighting_factors(:,:,i) = outflow_weighted(:,:,i).*...
            (outflow_weighted(:,:,i)>0)./outflow_gradients_sum;
    end


    %surface(z,outflow_weighting_factors(:,:,8))
    %figure
    %surface(z,mult_flow(:,:,3))
    % --------------------------
    % 4) The actual calculation of flow accumulation
    % temp_inflow_sum will store the amount of inflow in each iteration and
    % pass it on to total_flow_accumulation, which collects these inflows
    % until the final value is reached
    %shift_matrix = [1 0;1 -1;0 -1;-1 -1;-1 0;-1 1;0 1;1 1];
    temp_inflow_sum = outflow_gradients_sum*0+1;
    temp_outflow_weighting_factors = outflow_weighting_factors;
    total_flow_accumulation = outflow_gradients_sum*0+1;

    temp_inflow = gradients*0;

    %before:
    % Flow accumulation (terminates when there is no inflow):
    % fprintf(['processing topography for %i centroids...'],length(centroids.centroid_ID));
    while sum(temp_inflow_sum(:),'omitnan')>0
        for i = 1:8
            temp_inflow(:,:,i) = circshift(temp_inflow_sum,...
                [shift_matrix(i,1) shift_matrix(i,2)]);
        end
        % temp_inflow_sum collects the weighted current inflow from
        % contributing neighbouring cells (i.e. the ones with positive
        % gradients)
        temp_inflow_sum = sum(temp_inflow.*outflow_weighting_factors.*gradients>0,3);
        total_flow_accumulation = total_flow_accumulation + temp_inflow_sum;
    end




    % total_flow_accumulation(isnan(total_flow_accumulation))=0;fprintf(' done\n')
    % Calculate wetness index
    tmp_slope = slope + 0.1; %slope(slope==0) = min(min(slope(slope>0))); % we don't want -inf values for wet_index
    TWI = log((1+total_flow_accumulation)./tand(tmp_slope));
end


end


