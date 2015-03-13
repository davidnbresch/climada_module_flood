function centroids = centroids_fl_score_calc(centroids, check_plots, force_recalc)
% Calculate flood scores
% MODULE:
%   tbd
% NAME:
%	centroids_fl_score_calc
% PURPOSE:
%   Calculate flood score for given centroids
%   centroids_fl_score_calc applies a multiple-flow-direction method which, 
%   in contrast to the simple D8 method, allows the runoff to flow to 
%   multiple neighbouring cells. The distribution of the flow is determined 
%   based on the respective gradients between the central cell and its 
%   neighbouring cells. Also, the algorithm controls dispersion effects 
%   using a method suggested by 
%       Freeman, T.G. (1991): Calculating catchment area with 
%       divergent flow based on a regular grid; 
%       doi:10.1016/0098-3004(91)90048-I 
% CALLING SEQUENCE:
%   centroids = centroids_fl_score_calc(centroids, check_plots)
% EXAMPLE:
%   centroids = centroids_fl_score_calc(centroids_CHE,0)
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
%   centroids: centroids with an additional field 'flood_score', which 
%   assigns a number for flow accumulation to each centroid
% MODIFICATION HISTORY:
% Melanie Bieli, melanie.bieli@bluewin.ch, 20150226, initial
% Melanie Bieli, melanie.bieli@bluewin.ch, 20150311, added wetness index
%-

global climada_global

% check input arguments
if ~climada_init_vars; return; end
if ~exist('centroids','var') || isempty(centroids)
    cprintf([206 50 50]/255,['Error: Missing input ''centroids'' to '...
        'function centroids_fl_score_calc, can''t proceed.\n']);
    return;
end
if ~exist('check_plots', 'var') || isempty(check_plots) 
    check_plots = 0;    end
if ~exist('force_recalc','var') || isempty(force_recalc), 
    force_recalc = 0;  end


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
%

% We only calculate the basin IDs if the centroids do not come equipped
% with them
if ~isfield(centroids,'flood_score') || force_recalc
    
    % Get elevation data if centroids do not come equipped with them
    if ~isfield(centroids, 'elevation_m')
        cprintf([1,0.5,0],['\tWARNING: centroid elevation data missing, ' ...
            'using etopo.\n \tConsider running climada_read_srtm_DEM first.\n'])
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
            fprintf(['ERROR: no etopo_get function found. Please download the ' ...
                '<a href="https://github.com/davidnbresch/climada_module_etopo">'...
                'CLIMADA module etopo</a> from Github.\n'])
            return
        end
        TOPO_data=etopo_get(centroids_rect);
        if isempty(TOPO_data),return;end % error messages from etopo_get already
        centroids.elevation_m=interp2(TOPO_data.x,TOPO_data.y,TOPO_data.h,...
            centroids.lon,centroids.lat);
    end
    
    
    % Calculate mean difference between centroids in x and y direction
    % We assume a constant distance between degrees of latitude (lat), i.e. we
    % ignore the earth's slightly ellipsoid shape. The difference in
    % longitude (lon) is calculated based on:
    % length of 1 degree of lon = cosine(lat) * length of degree at equator
    [lon, lat] = meshgrid(unique(centroids.lon),unique(centroids.lat));
    
    dx = mean(diff(unique(centroids.lon)))*cos(mean(mean(lat))*pi/180)* 111.12 * 1000;
    dy = mean(diff(unique(centroids.lat)))* 111.12 * 1000;
    
    z       = griddata(centroids.lon,centroids.lat,centroids.elevation_m,lon,lat);
    c_ID    = griddata(centroids.lon,centroids.lat,centroids.centroid_ID,lon,lat);
    c_ID    = round(c_ID);
    
    % Calculate gradients in x and y direction in order to derive normal
    % vectors of the grid cells
    vx = zeros([size(z) 3]);
    vy = zeros([size(z) 3]);
    
    [dzdx, dzdy] = gradient(z, dx, dy);
    vx(:,:,3)=dzdx.*dx;     vx(:,:,1)=dx;
    vy(:,:,3)=dzdy.*dy;     vy(:,:,2)=dy;
    
    normal_vector=cross(vx,vy);
    area=sqrt(sum(normal_vector.^2,3));
    
    %Normalize Vector
    for ndx=1:3, normal_vector(:,:,ndx)=normal_vector(:,:,ndx)./area; end;
    
    % Calculate slope and aspect angle
    slope               =   acosd(normal_vector(:,:,3));
    aspect              =   atan2d(-1 .* normal_vector(:,:,2),normal_vector(:,:,1))+180;
    aspect(slope ==0)   =   nan;
    
    
    % Now we determine flow accumulation (which, in contrast to the local
    % parameters slope and aspect, can only be caluclated from the global
    % neighbourhood)
    % General concept: add cell inflows from higher adjacent cells, starting
    % from the the specified cell and working up to the watershed boundaries
    % The outflow of each cell is assigned to the neighbouring cell to which
    % the gradient is steepest
    
    % shift_matrix contains indices that will be used to access the
    % neighbouring cells
    shift_matrix = [0 -1;-1 -1;-1 0;1 -1;0 1;1 1;1 0;-1 1];
    [a, b]=size(z);
    % gradients are stored in a 3D-matrix - gradients(:,:,1) refers to the
    % eastern neighbour, gradients(:,:,2) to the southeastern neighbour etc.
    gradients = zeros(a,b,8);
    
    % For southeastern, southwestern, northwestern and northeastern neighbours,
    % the gradient is calculated by taking the elevation difference of the two
    % cells and dividing it by sqrt(2*grid_cell_size)
    % For the time being, we use dx as a mean value for grid cell size; for
    % NOAA's global 1 arc-minute ETOPO1 data, the horizontal resolution at the
    % equator would be approximately 1.85 kilometers):
    for c = 2:2:8
        average_grid_cell_size = (dx+dy)/2;
        gradients(:,:,c) = (circshift(z,[shift_matrix(c,1) shift_matrix(c,2)])...
            -z)/sqrt(2*average_grid_cell_size);
    end
    % For eastern, western, southern, and northern neighbours, the gradient is
    % simply calculated by taking the elevation difference of the two cells and
    % dividing it by the cell size
    
    for c = 1:2:7
        gradients(:,:,c) = ...
            (circshift(z,[shift_matrix(c,1) shift_matrix(c,2)])...
            -z)/dx;
    end
    gradients = atan(gradients)/pi*2;
    
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
    
    % --------------------------
    % 4) The actual calculation of flow accumulation
    % temp_inflow_sum will store the amount of inflow in each iteration and
    % pass it on to total_flow_accumulation, which collects these inflows
    % until the final value is reached
    temp_inflow_sum = outflow_gradients_sum;
    total_flow_accumulation = outflow_gradients_sum;
    
    temp_inflow = gradients*0;
    
    % Flow accumulation (terminates when there is no inflow):
    fprintf('Calculating flood scores for %i centroids...\n',...
        length(centroids.centroid_ID));
    while sum(temp_inflow_sum(:))>0
        % iterate over inflow from all directions and store it in temp_inflow
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
    total_flow_accumulation(isnan(total_flow_accumulation))=0;
    
    % Calculate wetness index
    slope(slope==0) = min(min(slope(slope>0))); % we don't want inf values for wet_index
    wet_index = log((1+total_flow_accumulation)./tand(slope));
    
    % ---------------------------
    % Add field flood_score to the centroids struct, i.e. loop over all
    % centroid IDs (which do not in all cases start at 1!) and assign the
    % respective flow accumulation numbers (the so-called flood scores)
    fprintf('Assigning flood scores to centroids...\n')
    for centroid_i = centroids.centroid_ID(1):centroids.centroid_ID(end)
        centroids_indices = find(c_ID == centroid_i);
        centroids.flood_score(centroid_i) = ...
            mean(total_flow_accumulation(centroids_indices));
        centroids.wetness_index(centroid_i) = ...
            mean(wet_index(centroids_indices));
    end
    centroids.flood_score(isnan(centroids.flood_score))=0;
    centroids.wetness_index(isnan(centroids.wetness_index))=0;
    % set flood scores and wetness indices to 0 for centroids in the ocean
    centroids.flood_score(centroids.onLand==0)=0;
    centroids.wetness_index(centroids.onLand==0)=0;
    
    cprintf([23 158 58]/255,['Successfully completed calculation of '...
        'flood scores.\n']);
else
    % centroids already have a field 'flood_scores'
    cprintf([23 158 58]/255,'Skipped - centroids already have flood scores.\n')
end % if isfield(centroids,'flood_score')
    

% -----------------------------
% If requested, generate pseudocolor plots of elevation, slope, aspect 
% angle, and flow accumulation 
if check_plots
    
    % plot elevation data
    figure('Name','Elevation','Color',[1 1 1]);
    dem(unique(lon),unique(lat),z)
    
    % plot slope
    figure('Name','Slope','Color',[1 1 1]);    
    h = pcolor(slope);
    colormap(jet), colorbar
    set(h,'LineStyle','none')
    axis equal
    title('Slope [degrees]')
    [r, c] = size(slope);
    axis([1 c 1 r])
    set(gca,'TickDir','out')
  
    
    % plot aspect angle
    figure('Name','Aspect','Color',[1 1 1]);
    h=pcolor(aspect);
    colormap(hsv),colorbar
    set(h,'Linestyle','none')
    axis equal
    title('Aspect')
    axis([1 c 1 r])
    set(gca,'TickDir','out')
    
    % Plot flow accumulation
    figure('Name','Flood scores','Color',[1 1 1]);
    h = pcolor(log(1+total_flow_accumulation));
    %h = pcolor(total_flow_accumulation);
    colormap(flipud(jet)), colorbar
    set(h,'LineStyle','none')
    axis equal
    title('Flood scores')
    [r, c] = size(total_flow_accumulation);
    axis([1 c 1 r])
    set(gca,'TickDir','out')
    
    % Plot wetness index
    figure('Name','Wetness index','Color',[1 1 1]);
    h = pcolor(wet_index);
    colormap(flipud(jet)), colorbar
    set(h,'LineStyle','none')
    axis equal
    title('Wetness index')
    [r, c] = size(wet_index);
    axis([1 c 1 r])
    set(gca,'TickDir','out')

    
%   % plot flow direction
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

end % if check_plots

end


