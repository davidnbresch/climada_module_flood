function centroids = climada_centroids_TWI_calc_v2(centroids,centroids_set_file,check_plots)
% Calculate flood scores and topographic wetness indices
% MODULE:
%   flood
% NAME:
%	climada_centroids_TWI_calc
% PURPOSE:
%   Calculate flood scores and topographic wetness indices for given
%   centroids, that are on a regular grid.
%   climada_centroids_TWI_calc applies a multiple-flow-direction method which,
%   in contrast to the simple D8 method, allows the runoff to flow to
%   multiple neighbouring cells. The distribution of the flow is determined
%   based on the respective gradients between the central cell and its
%   neighbouring cells. Also, the algorithm controls dispersion effects
%   using a method suggested by
%       Freeman, T.G. (1991): Calculating catchment area with
%       divergent flow based on a regular grid;
%       doi:10.1016/0098-3004(91)90048-I
%       https://ac.els-cdn.com/009830049190048I/1-s2.0-009830049190048I-main.pdf?_tid=866d08a6-c54f-11e7-b499-00000aab0f01&acdnat=1510233280_58e175c46862affef6785a80491025fe
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
%   centroids: Climada centroids struct; the following fields are required:
%         .lat:           Latitude
%         .lon:           Longitude
%         .centroid_ID:   centroid ID
%         .admin0_name    Country name
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
% Melanie Bieli, melanie.bieli@bluewin.ch, 20150226, initial
% Melanie Bieli, melanie.bieli@bluewin.ch, 20150311, added wetness index
% Gilles Stassen, gillesstassen@hotmail.com, 20150407, clean up
% Lea Mueller, muellele@gmail.com, 20150720, bugfix, quick+dirty workaround
%              to create meshgrid and allocate FL_score, @Gilles: please check and correct
% Jacob Anz, 280715, fixed shift_matrix
% Lea Mueller, muellele@gmail.com, 20150925, add process management/waitbar
% Lea Mueller, muellele@gmail.com, 20151105, improve output documentation
% Lea Mueller, muellele@gmail.com, 20151125, rename to climada_centroids_TWI_calc from centroids_TWI
% David N. Bresch, david.bresch@gmail.com, 20170629, double() of single for griddata
% Thomas R�lli, thomasroelli@gmail.com, 20180222, remove interpolation, and
%   take values from centroids.elevation_m 
% Thomas R�lli, thomasroelli@gmail.com, 20180223, init version 2

global climada_global

% check input arguments
if ~climada_init_vars; return; end
if ~exist('centroids',  'var') || isempty(centroids),   climada_centroids_load; end
if ~exist('check_plots','var') || isempty(check_plots), check_plots = 0;        end

% PARAMETERS
%
% To allow for smooth interpolation, we will enlarge the rectangle defined
% by the min and max latitude and longitude coordinates of the centroids by
% a certain amount (in degrees) at all four sides
add_degree = 1;
%
% Weighting factor to be applied in the calculation of flow accumulation
% according to Freeman (1991)
weighting_factor = 1.1;

% Get elevation data if centroids do not come equipped with them
if ~isfield(centroids, 'elevation_m')
    cprintf([1,0.5,0],['WARNING: centroid elevation data missing, ' ...
        'using etopo.\n \tConsider running climada_90m_DEM for ultra high resolution topography.\n'])
    % prep the region we need (rectangular region encompassing the hazard
    % centroids)
    centroids_rect=[min(centroids.lon) max(centroids.lon) ...
        min(centroids.lat) max(centroids.lat)];
    % enlarge the rectangle
    centroids_rect=[centroids_rect(1)-add_degree centroids_rect(2)+add_degree ...
        centroids_rect(3)-add_degree centroids_rect(4)+add_degree];
    
    % cut the elevation data out of the global topography dataset
    if ~exist('etopo_get','file')
        % safety to inform the user in case he misses the ETOPO module
        fprintf(['ERROR: no etopo_get function found. Please download ' ...
            '<a href="https://github.com/davidnbresch/climada_module_elevation_models">'...
            'elevation_models</a> from Github.\n'])
        return
    end
    TOPO_data=etopo_get(centroids_rect);
    if isempty(TOPO_data),return;end % error messages from etopo_get already
    centroids.elevation_m=interp2(TOPO_data.x,TOPO_data.y,TOPO_data.h,...
        centroids.lon,centroids.lat);
end

centroids.lon=double(centroids.lon); % to double, as required by griddata below
centroids.lat=double(centroids.lat);
centroids.elevation_m=double(centroids.elevation_m);
    
% Calculate mean difference between centroids in x and y direction
% We assume a constant distance between degrees of latitude (lat), i.e. we
% ignore the earth's slightly ellipsoid shape. The difference in
% longitude (lon) is calculated based on:
% length of 1 degree of lon = cosine(lat) * length of degree at equator
% lon_singleton = [min(centroids.lon):min(diff(unique(centroids.lon))):max(centroids.lon)];
% lat_singleton = [min(centroids.lat):min(diff(unique(centroids.lat))):max(centroids.lat)];
% 20180222 remove interpolation (because of slight change of coordinates
% and elevation
% removed:
%   lon_singleton = [min(centroids.lon):diff(centroids.lat(1:2))*factor_f:max(centroids.lon)];
%   lat_singleton = [min(centroids.lat):diff(centroids.lat(1:2))*factor_f:max(centroids.lat)];
%   [lon, lat] = meshgrid(lon_singleton,lat_singleton);
%   lon=double(lon);lat=double(lat); % to double, as required by griddata below
n_lon = numel(unique(centroids.lon));
n_lat = numel(unique(centroids.lat));
lon = reshape(centroids.lon,n_lat,n_lon);
lat = reshape(centroids.lat,n_lat,n_lon);

lon_singleton = lon(1,:);
lat_singleton = lat(:,1)';


% assume cos(lat) doesn't vary much within small study region and take cos
% of mean latitude
x = lon_singleton .* (cos(mean(mean(lat))*pi/180)* 111.12 * 1000); 

y = lat_singleton .* (111.12 * 1000);

factor_f = 1;

dx_ = min(diff(unique(centroids.lon))) * (cos(mean(mean(lat))*pi/180)* 111.12 * 1000);
dy_ = min(diff(unique(centroids.lat))) * (111.12 * 1000);
dx = diff(centroids.lat(1:2))*factor_f * (cos(mean(mean(lat))*pi/180)* 111.12 * 1000);
dy = diff(centroids.lat(1:2))*factor_f * (111.12 * 1000);

% dx = mean(diff(unique(centroids.lon)))*cos(mean(mean(lat))*pi/180)* 111.12 * 1000;
% dy = mean(diff(unique(centroids.lat)))* 111.12 * 1000;

% 20180222: removed
%    z       = griddata(centroids.lon,centroids.lat,centroids.elevation_m,lon,lat, 'linear');
%    c_ID    = griddata(centroids.lon,centroids.lat,centroids.centroid_ID,lon,lat, 'nearest');

z = reshape(centroids.elevation_m,n_lat,n_lon);
c_ID = reshape(centroids.centroid_ID,n_lat,n_lon);

if dx_ ~= dx || dy_ ~= dy || numel(z) > numel(centroids.centroid_ID)
    cprintf([1 0.5 0],'WARNING: centroids not defined by uniform rectangular grid - code will continue, but may encounter issues. Please check fields:\n')
    cprintf([1 0.5 0],'\t\t consider generating centroids on a uniform grid before running climada_centroids_TWI_calc\n')
    cprintf([1 0.5 0],'\t\t see climada_centroids_generate\n')
end

% tmp     = [num2str(centroids.lon') num2str(centroids.lon')]; 
% Calculate gradients in x and y direction in order to derive normal
% vectors of the grid cells
vx = zeros([size(z) 3]);
vy = zeros([size(z) 3]);

[dzdx, dzdy] = gradient(z, x, y);
vx(:,:,3)=dzdx.*dx;     vx(:,:,1)=dx;
vy(:,:,3)=dzdy.*dy;     vy(:,:,2)=dy;

normal_vector   =   cross(vx,vy);
area            =   sqrt(sum(normal_vector.^2,3));

%Normalize Vector
for ndx=1:3, normal_vector(:,:,ndx)=normal_vector(:,:,ndx)./area; end;

% Calculate slope and aspect angle
slope               =   acosd(normal_vector(:,:,3));
aspect              =   atan2d(-1 .* normal_vector(:,:,2),normal_vector(:,:,1))+180;
aspect(slope ==0)   =   nan;

%%%% the following works well only for a regular grid of centroids
% find sink centroid
flow_dir            =   round(aspect./45);  % discretise aspect to 1 of 8 possible directions (neighbouring cells)
flow_dir(flow_dir == 0)   = 8;              % direction 8 is equivalent to 0;
flow_dir(isnan(aspect))   = 0;              % redefine zeroth direction as no outward flow

sink_cID                =   c_ID; % init

% direction unit vectors [e.g. flow_dir = 3 => vector = (di(3),dj(3))]
di = [-1  0  1  1  1  0 -1 -1]; 
dj = fliplr([ 1  1  1  0 -1 -1 -1  0]);

[n_j,n_i] = size(c_ID);
i_ndx = [1:n_i]; j_ndx = [1:n_j];

for dir_i = 1:8
    % construct index reordering vector for i direction
    tmp_i_ndx                   = i_ndx + di(dir_i);
    % sink of boundary cell is itself, use expression after '%' for periodic
    % boundary conditions
    tmp_i_ndx(tmp_i_ndx <= 0)   = 1;    % n_i - tmp_i_ndx(tmp_i_ndx <= 0);
    tmp_i_ndx(tmp_i_ndx > n_i)  = n_i;  % tmp_i_ndx(tmp_i_ndx > n_i) - n_i;
    
    % construct index reordering vector for i direction
    tmp_j_ndx                   = j_ndx - dj(dir_i);
    % sink of boundary cell is itself, use expression after '%' for periodic
    % boundary conditions
    tmp_j_ndx(tmp_j_ndx <= 0)   = 1;    % n_j - tmp_j_ndx(tmp_j_ndx <= 0);
    tmp_j_ndx(tmp_j_ndx > n_j)  = n_j;  % tmp_j_ndx(tmp_j_ndx > n_j) - n_j;
    
    tmp_c_ID    = c_ID(tmp_j_ndx,tmp_i_ndx);
    sink_cID(flow_dir == dir_i)  = tmp_c_ID(flow_dir == dir_i);
end
sink_cID(flow_dir == 0)     =   NaN;
sink_cID(sink_cID == c_ID)  =   NaN;
clear tmp_c_ID tmp_i_ndx tmp_j_ndx


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
shift_matrix = [0 -1;-1 -1;-1 0;-1 1;0 1;1 1;1 0;1 -1]; %Thomas R�lli 20180223 [c(8,:) and c(4,:) was the same before]
[a, b]=size(z);
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
    gradients(:,:,c) = (circshift(z,[shift_matrix(c,1) shift_matrix(c,2)])...
        -z)/cell_diag;
end
% For eastern, western, southern, and northern neighbours, the gradient is
% simply calculated by taking the elevation difference of the two cells and
% dividing it by the cell size

for c = 1:2:7
    if c == 1 || c == 5, distance=dx; else distance=dy; end
    gradients(:,:,c) = (circshift(z,[shift_matrix(c,1) shift_matrix(c,2)])-z)/distance;
end

%lat and z need to be flipped updown (requested by
%climada_centroids_gradients) afterwards reverse it to have same structure
%again.
% [gradients,~,~] = climada_centroids_gradients(flipud(lon),flipud(lat),flipud(z));
% gradients = flipud(atand(gradients));

% set gradients which are directed out of the study region to zero; at
% boarder cells

gradients(numel(lat(:,1)),:,[2 3 4]) = 0; %upper boarder
gradients(1,:,[6 7 8]) = 0; %lower boarder
gradients(:,1,[4 5 6]) = 0; %left boarder
gradients(:,numel(lon(1,:)),[1 2 8]) = 0; %right boarder
gradients = atand(gradients);%/pi*2;

[gradients2,~,~] = climada_centroids_gradients(flipud(lon),flipud(lat),flipud(z));
gradients2 = flipud(atand(gradients2));

%----------------------
% Prep for the calculation of flow accumulation (resulting in flood
% scores), where the outflow of each cell will be distributed among its
% neighbouring cells based on the steepness of the respective gradients
% 
% 1) Take the absolute value of negative gradients (i.e., the ones that
% cause outflow) and raise them to the power of the selected weighting
% factor.
% Flow of positive gradients is set to zero.
%gradients = gradients*-1;
%gradients(gradients < 0) = 0; 
load('C:\Users\Simon R�lli\Desktop\data\mult_flow.mat','mult_flow');
%load('C:\Users\Simon R�lli\Desktop\data\gradients2.mat','gradients2');
load('C:\Users\Simon R�lli\Desktop\data\wet_index_wrong.mat','wet_index_wrong')
% for c=1:8
%     gradients2(:,:,c) = atand(flipud(gradients2(:,:,c)));
%     mult_flow(:,:,c) = flipud(mult_flow(:,:,c));
% end

%weighting_factor = 4;
outflow_weighted = (gradients.*(-1*gradients<0)).^weighting_factor;
%outflow_weighted = (gradients.*(gradients<0)*-1).^weighting_factor;
outflow_weighted2 = (gradients2.*(gradients2<0)*-1).^weighting_factor;
% outflow_weighted = (abs(gradients)).^weighting_factor;

% 2) Sum up gradients such that for each grid cell, outflow_gradients_sum
% contains the sum of the gradients to the grid cell's 8 neighbouring cells
outflow_gradients_sum = sum(outflow_weighted,3);
outflow_gradients_sum2 = sum(outflow_weighted2,3);
outflow_gradients_sum(outflow_gradients_sum==0) = 1; %prevent division by 0
outflow_gradients_sum2(outflow_gradients_sum2==0) = 1; %prevent division by 0

% 3) Normalize the flow such that the neighboring cell gradients of each
% direction (east, southeast, south, etc.) add up to 1. These fractions
% will be used to distribute the flow to all neighbours of a cell
for i = 1:8
    outflow_weighting_factors(:,:,i) = outflow_weighted(:,:,i).*...
        (outflow_weighted(:,:,i)>0)./outflow_gradients_sum;
    outflow_weighting_factors2(:,:,i) = outflow_weighted2(:,:,i).*...
        (outflow_weighted2(:,:,i)>0)./outflow_gradients_sum2;
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
count=0;

% figure('units','normalized','outerposition',[0 0 1 1])
% s = surf(lon,lat,z,total_flow_accumulation);
% colorbar
%caxis([10000 11000])

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
%     s.CData = total_flow_accumulation;
%     pause(0.05)
end


% temp_inflow_sum = outflow_gradients_sum2*0+1;
% temp_outflow_weighting_factors = outflow_weighting_factors2;
% total_flow_accumulation2 = outflow_gradients_sum2*0+1;
% figure('units','normalized','outerposition',[0 0 1 1])
% s = surf(lon,lat,z,total_flow_accumulation2);
% colorbar
% %lim = caxis
% %caxis([0 1200])
% 
% temp_outflow = gradients*0;
% count2=0;

% shift_matrix = [1 0;1 1;0 1;-1 1;-1 0;-1 -1;0 -1;1 -1];
% %before:
% % Flow accumulation (terminates when there is no inflow):
% % fprintf(['processing topography for %i centroids...'],length(centroids.centroid_ID));
% ud = [5 4 3 2 1 8 7 6];
% outflow_weighting_factors2_ud = outflow_weighting_factors2;
% for c=8
%     outflow_weighting_factors2_ud(:,:,c)=outflow_weighting_factors2(:,:,ud(c));
% end
% 
% while sum(temp_inflow_sum(:),'omitnan')>0
%     % iterate over inflow from all directions and store it in temp_inflow
%     for i = 1:8
%         %temp_inflow(:,:,i) = circshift(temp_inflow_sum,...
%             %[shift_matrix(i,1) shift_matrix(i,2)]).*outflow_weighting_factors2(:,:,i);
%         %temp_outflow(:,:,i) = temp_inflow_sum.*(outflow_weighting_factors2(:,:,i)>0);
%         temp_inflow(:,:,i) = circshift(temp_inflow_sum,[shift_matrix(i,1) shift_matrix(i,2)])...
%             .*circshift(outflow_weighting_factors2_ud(:,:,i),[shift_matrix(i,1) shift_matrix(i,2)]).*(gradients(:,:,i)>0);
%         %temp_outflow(:,:,i) = temp_inflow_sum.*(outflow_weighting_factors2(:,:,i));
%         %temp_inflow(:,:,i) = circshift(temp_outflow(:,:,i),[shift_matrix(i,1) shift_matrix(i,2)]);
% %         figure
% %         surface(z,temp_outflow(:,:,i))
% %         figure
% %         surface(z,temp_inflow(:,:,i))
%            
%     end
%     % temp_inflow_sum collects the weighted current inflow from
%     % contributing neighbouring cells (i.e. the ones with positive
%     % gradients)
%     %sum(temp_inflow_sum)
%     count2=count2+1
%     temp_inflow_sum = sum(temp_inflow,3);
%     sum(sum(temp_inflow_sum))
%     total_flow_accumulation2 = total_flow_accumulation2 + temp_inflow_sum;
%     s.CData = total_flow_accumulation2;
%     
%     pause(0.075)
% end

% total_flow_accumulation(isnan(total_flow_accumulation))=0;
fprintf(' done\n')
% Calculate wetness index
tmp_slope = slope + 0.1; %slope(slope==0) = min(min(slope(slope>0))); % we don't want -inf values for wet_index
wet_index = log((1+total_flow_accumulation)./tand(tmp_slope));
%wet_index2 = log((1+total_flow_accumulation2)./tand(tmp_slope));
% figure
% surface(z,wet_index)
% figure
%surface(z,wet_index2)

% ---------------------------
% Add field flood_score to the centroids struct, i.e. loop over all
% centroid IDs (which do not in all cases start at 1!) and assign the
% respective flow accumulation numbers (the so-called flood scores)
% fprintf('assigning derived topographic properties to centroids...')

% process managemagent
msgstr   = sprintf('Assigning topographic properties to %i centroids ... ',centroids.centroid_ID(end));
mod_step = 10; % first time estimate after 10 assets, then every 100
if climada_global.waitbar
    fprintf('%s (updating waitbar with estimation of time remaining every 100th centroid)\n',msgstr);
    h        = waitbar(0,msgstr);
    set(h,'Name','Assigning TWI');
else
    fprintf('%s (waitbar suppressed)\n',msgstr);
    format_str='%s';
end

% init
centroids.FL_score   = zeros(size(centroids.lon));
centroids.TWI        = zeros(size(centroids.lon));
centroids.slope_deg  = zeros(size(centroids.lon));
centroids.area_m2    = zeros(size(centroids.lon));
centroids.aspect_deg = zeros(size(centroids.lon));
centroids.sink_ID    = zeros(size(centroids.lon));

for centroid_i = centroids.centroid_ID(1):centroids.centroid_ID(end)
    centroids_ndx                   = c_ID == centroid_i;    
    if sum(centroids_ndx(:))==1 
        
        centroids.area_m2(centroid_i) 	= area      (centroids_ndx);
        centroids.aspect_deg(centroid_i)= aspect    (centroids_ndx);
        centroids.FL_score(centroid_i)  = total_flow_accumulation(centroids_ndx);
        centroids.sink_ID(centroid_i)   = sink_cID  (centroids_ndx);
        centroids.slope_deg(centroid_i) = slope     (centroids_ndx);
        centroids.TWI(centroid_i)       = wet_index (centroids_ndx);
    else
        centroids.area_m2(centroid_i) 	= NaN;
        centroids.aspect_deg(centroid_i)= NaN;
        centroids.FL_score(centroid_i)  = NaN;
        centroids.sink_ID(centroid_i)   = NaN; %mode(mode(sink_cID(centroids_ndx)));
        centroids.slope_deg(centroid_i) = NaN;
        centroids.TWI(centroid_i)       = NaN;
    end
    % the progress management
    if mod(centroid_i,mod_step)==0
        mod_step = 100;
        msgstr = sprintf('%i/%i centroids',centroid_i,centroids.centroid_ID(end));
        if climada_global.waitbar
            waitbar(centroid_i/centroids.centroid_ID(end),h,msgstr); % update waitbar
        else
            fprintf(format_str,msgstr); % write progress to stdout
            format_str=[repmat('\b',1,length(msgstr)) '%s']; % back to begin of line
        end
    end
end

% centroids.FL_score(isnan(centroids.FL_score))=0;
% centroids.TWI(isnan(centroids.TWI))=0;
% set flood scores and wetness indices to 0 for centroids in the ocean
% and in the buffer zone
centroids.FL_score(centroids.onLand==0)=NaN; %0;
centroids.FL_score(centroids.onLand==2) = NaN;%0;
centroids.TWI(centroids.onLand==0)=NaN; %0;
centroids.TWI(centroids.TWI<0)=NaN; %0;
%centroids.TWI(centroids.onLand==2)) = 0;

fprintf(' done\n');

% -----------------------------
% If requested, generate pseudocolor plots of elevation, slope, aspect
% angle, and flow accumulation
if check_plots
    
%     % plot elevation data
%     figure('Name','Elevation','Color',[1 1 1]);
%     climada_DEM_plot(unique(lon),unique(lat),z)
%     
%     % plot slope
%     figure('Name','Slope','Color',[1 1 1]);
%     h = pcolor(slope);
%     colormap(jet), colorbar
%     set(h,'LineStyle','none')
%     axis equal
%     title('Slope [degrees]')
%     [r, c] = size(slope);
%     axis([1 c 1 r])
%     set(gca,'TickDir','out')
%     
%     % plot aspect angle
%     figure('Name','Aspect','Color',[1 1 1]);
%     h=pcolor(aspect);
%     colormap(hsv),colorbar
%     set(h,'Linestyle','none')
%     axis equal
%     title('Aspect')
%     axis([1 c 1 r])
%     set(gca,'TickDir','out')
%     
%     % Plot flow accumulation
%     figure('Name','Flood scores','Color',[1 1 1]);
%     h = pcolor(log(1+total_flow_accumulation));
%     %h = pcolor(total_flow_accumulation);
%     colormap(flipud(jet)), colorbar
%     set(h,'LineStyle','none')
%     axis equal
%     title('Flood scores')
%     [r, c] = size(total_flow_accumulation);
%     axis([1 c 1 r])
%     set(gca,'TickDir','out')
    
    % Plot wetness index
    figure('Name','Topographic wetness index','Color',[1 1 1]);
    h = pcolor(wet_index);
    colormap(flipud(jet)), colorbar
    set(h,'LineStyle','none')
    axis equal
    title('Topographic wetness index')
    [r, c] = size(wet_index);
    axis([1 c 1 r])
    set(gca,'TickDir','out')
end

% plot flow direction
%   figure('Name','Rainfall runoff direction', 'Color', [1 1 1])
%    title('Runoff flow direction')
%   h = contourf(z, linspace(min(min(z)),max(max(z)),15));
%   hold on
%   gradients(isinf(gradients)) = nan;
%   quiver(x,y, dzdx.*-1, dzdy.*-1);
%   %quiver(x, y, gradients(:,:,1)+gradients(:,:,5), ...
%     %gradients(:,:,3)+gradients(:,:,7));
%   %[size_dim1, size_dim2, ~] = size(gradients);
%   %axis([1 size_dim1 1 size_dim2])
%   %set(gca,'TickDir','out')
%   hold off

%   % plot runoff direction
%   figure('Name','Rainfall runoff direction', 'Color', [1 1 1])
%   title('Runoff flow direction')
%   hold on
%   contourf(x,y,z)
%   quiver(x, y, dzdx.*-1, dzdy.*-1)