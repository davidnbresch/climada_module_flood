function [rel_prob,intensity] = climada_ls_propagation(source_area,mult_flow,... 
hor_dist,ver_dist,v_max,phi,delta_i,perWt,en_decay,cal_int,cell_area,stdv_wiggle)

% Computes the flow path according to the multiple flow algorithm
% (according to Holmgren 1994). The flow distance is taken into account by
% a simplified friction model (see Horton et al. 2013).
% MODULE:
%   flood
% NAME:
%   climada_ls_propagation
% PURPOSE:
%   
% CALLING SEQUENCE:
%   climada_ls_propagation(source_area,mult_flow,hor_dist,ver_dist,v_max,phi,delta_i,perWt)
% EXAMPLE:
%   
% INPUTS: 
%   source_area: matrix with information about the source areas (0/1).
%                The quantity of the 3rd dimension of the matrix
%                corresponds to the number of events
%   mult_flow: matrix (same nxm dimensions as source_areas) with outflow 
%              proportion of each cell to its neighbours.
%              Dimension: nxmx8
%   hor_dist: horizontal distance of each cell to its 8 neighbours.
%             Dimension:nxmx8. Needed to calculate friction
%   ver_dist: vertical distance of each cell to its 8 neighbours.
%             Dimension:nxmx8. Needed to calculate potential energy
% OPTIONAL INPUT PARAMETERS:
%   v_max:   describes maximal possible velocity of slide. If velocity is 
%            exceeded v_max is taken. Should keep energy amounts within reasonal values
%            and therefore prevent improbalbe runout distances
%   phi:     angle of reach, angle of the line connnecting the source area to
%            the most distant point reached by the slide, along its path.
%            Factor controlls maximum possible runout distance
%   delta_i: Minimum threshold which prevent propagation when deceeded. Small values
%            are removed and spreaded to the other remaining, lower
%            situated neighbouring cells 
%   perWt:   row vector with 8 elements. gives weight to each direction.
%            The persistence function aims at
%            reproducing the behaviour of inertia --> weights the flow
%            direction according to change of direction to neighbour.
%            First element represent weight to neighbour in same direction
%            as flow (0 degree), second element weights right neighbour 45
%            degrees (angle between previous direction and direction from
%            the central cell to corresponding neighbour) .. last element
%            represents weight of neighbour to the left of flow direction
%            (45 degree anticlockwise). 
%   en_decay:fraction of energy which shall be removed from total energy
%            after a cell was propagated to the next cell -> pragmatic 
%            approach to reduce the number of too long slides
%   cal_int:     (0/1) set to 1 if intenisty should be calculated (in m^2),
%            Intensity is defined as covered area (similar to upstream
%            contributing area). If one cell is passed more than once
%            , the maximum intensity is saved. If set 0 (default): a zeromatrix
%            is returned.
%   cell_area: need to be provided if cal_int set to 1. derived from
%              climada_centroids_area
%   stdv_wiggle: Matrix(nlat x nlon) of wiggle factor for each source cell
%                as a product of the stdv of the slope degree (provided by
%                sqrt(climada_centroids_demVariance)) and a random normal
%                distr. number (provided by abs(normrnd(0,1))). If stdv_wiggle
%                is empty it is set to 0 (default). If not empty: slope at
%                source cell is wiggled by slope(x,y) + stdv_wiggle (each
%                source cell).
%                A positive wiggle increases the downward slope gradient by
%                the same amount.
%                positive slopes to neighbours (upward slope) are not affected by wiggle
%                --> to prevent upward slope --> potential energy is
%                calculated by taking slope which is not wiggled
%                Method is implemented to prevent too many slopes which are
%                not triggered because of too flat area. Therefore the
%                standarddeviation should be multiplied with a positive
%                normal distibuted random number (abs(normrnd(0,1)))
% OUTPUTS:
%   rel_prob: matrix (lon lat) with the resulting relative probability after a
%                  single landslides is propagated downstream starting from source areas
%   intensity:   matrix (lon lat) which saves the intensity [m^2] as
%                  covered area (similar to upstream contributing area)
%  
% MODIFICATION HISTORY:
% Thomas Rölli, thomasroelli@gmail.com, 20180219, init
% Thomas Rölli, thomasroelli@gmail.com, 20180227, changed shift_matrix
% Thomas Rölli, thomasroelli@gmail.com, 20180305, implement v_max of slide
% Thomas Rölli, thomasroelli@gmail.com, 20180317, init v2
% Thomas Rölli, thomasroelli@gmail.com, 20180319, changed structure, add
%  persistence function, add intensity threshold (delta_i)
% Thomas Rölli, thomasroelli@gmail.com, 20180329, sum of intensity while
%  propagated
% Thomas Rölli, thomasroelli@gmail.com, 20180404, renamed from
%  climada_ls_spread_v2 to climada_ls_propagation
% Thomas Rölli, thomasroelli@gmail.com, 20180406, calculate distance to
%  source
% Thomas Rölli, thomasroelli@gmail.com, 20180529, use of indices for active
%  cells instead of matrices --> faster
% Thomas Rölli, thomasroelli@gmail.com, 20180608, implementation of
%  intensity definition (as covered area)
% Thomas Rölli, thomasroelli@gmail.com, 20180629, implementation of
%  pragmatic reduction of energy
% Thomas Rölli, thomasroelli@gmail.com, 20180718, include wiggle of source
%  cell according to standarddeviation of slope 

global climada_global
if ~climada_init_vars, return; end

% check arguments
if ~exist('source_area', 'var'), source_area = []; end
if ~exist('mult_flow', 'var'), mult_flow = []; end
if ~exist('hor_dist', 'var'), hor_dist = []; end
if ~exist('ver_dist', 'var'), ver_dist = []; end
if ~exist('v_max', 'var'), v_max = []; end
if ~exist('phi', 'var'), phi = []; end
if ~exist('delta_i', 'var'), delta_i = []; end
if ~exist('perWt', 'var'), perWt = []; end
if ~exist('en_decay', 'var'), perWt = []; end
if ~exist('cal_int', 'var'), cal_int = []; end
if ~exist('cell_area', 'var'), cell_area = []; end
if ~exist('stdv_wiggle', 'var'), stdv_wiggle = []; end

% PARAMETERS 
if isempty(source_area); return; end
if isempty(mult_flow); return; end
if isempty(hor_dist); return; end
if isempty(ver_dist); return; end
if isempty(v_max); v_max = 4; end
if isempty(phi); phi = 27; end %empirical minimum travel angle, used for friction-calculation
if isempty(delta_i); delta_i = 0.0003; end
if isempty(perWt); perWt = [1 0.8 0.4 0 0 0 0.4 0.8]; end
if isempty(en_decay); en_decay = 0; end
if isempty(cal_int); cal_int = 0; end
if isempty(cell_area) && cal_int==1; return; end
if isempty(stdv_wiggle); stdv_wiggle = zeros(size(source_area)); end

%if numel(perWt) ~= 8; return; end

%for calculations
g = 9.81; %acceleration of gravity 

n_lat = numel(source_area(:,1));
n_lon = numel(source_area(1,:));

%assessing flow path; starting from cells which are equal 1 in
spread = double(source_area);
active_cells = source_area; %save active cells for next iteration
temp_active_cells = active_cells*0;
energy = double(zeros(n_lat,n_lon)); %save energy of cells while propagating
eFric = g*hor_dist*tand(phi); %loss of energy by friction, for whole field in each direction
ePot = g*ver_dist*(-1); %gain of potential enenergy, for whole field in each direction
shift_matrix = [-1 0;-1 -1;0 -1;1 -1;1 0;1 1;0 1;-1 1]*-1;
active_cells_idx = []; %saves the indices of all active cells 


%direction-matrix --> saves direction of flow (1-8)in each iteration
direction = source_area*0; %will save flow-direction in each iteration

%total intensity
rel_prob = double(source_area);

%distance to source
intensity = zeros(size(source_area));
if cal_int
    %set intensity of source cell to area of source cell
    intensity(source_area==1) = cell_area(source_area==1);
end

%%%%%for animation during propagation can be removed%%%%%
% load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_srtm1_centroids.mat')
% n_lon = numel(unique(centroids.lon));
% n_lat = numel(unique(centroids.lat));
% lon = reshape(centroids.lon,n_lat,n_lon);
% lat = reshape(centroids.lat,n_lat,n_lon);
% 
% picnumber = 1;
% path = 'C:\Users\Simon Rölli\Desktop\climada\climada\animation';
% fullpath = [path filesep 'pic' num2str(picnumber) '.tif'];
% figure('units','normalized','outerposition',[0 0 1 1])
% 
% prob_plot = rel_prob;
% prob_plot(prob_plot==0) = nan;
% s = surface(lon,lat,rel_prob);
% s.CData = log(prob_plot);
% K = GRIDobj(lon,lat,log(prob_plot));
% GRIDobj2geotiff(K,fullpath);
% colorbar
% caxis([-4 0])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k=0;
while sum(sum(active_cells))>0 %iteration through number of runs  
temp_active_cells = temp_active_cells*0;
active_cells_idx = find(active_cells == 1);
for k=1:numel(active_cells_idx)
    %get j and i (should be switches --> i for row, j for column)
    [j,i] = ind2sub(size(source_area),active_cells_idx(k));
    %if active_cells(j,i)
        %extract information out of mult_flow
        temp_mult_flow = mult_flow(j,i,:);
        temp_mult_flow = temp_mult_flow(:)';
        %check if flow just started (source area) --> if yes no
        %influence of persistence function
        if direction(j,i) ~= 0 
            shifted_perWt = circshift(perWt,direction(j,i)-1);
        else
            shifted_perWt = ones(1,8);
        end
        %suceptibility = multiple flow * persistence function
        wgt_suscept = temp_mult_flow.*shifted_perWt;
        %scale such that sum equal 1 --> outflow proportion in each
        %direction
        wgt_suscept = wgt_suscept./sum(wgt_suscept);

        %calculate energy to all its neigbours
        if direction(j,i) == 0 %source cell --> do not remove energy by en_decay but wiggle start slope by stdv_wiggle
            %calculate poential energy if slope is wiggled (for all 8
            %neighbours)
            slope = atand(ver_dist(j,i,:)./hor_dist(j,i,:)); %slope if not wiggled
            wiggled_slope = slope*(-1) + stdv_wiggle(j,i); %wiggled slope: *(-1) to have downslope positive --> positive wiggle lead to steeper downslope
            %pot energy to neigbour with wiggled slope --> hor_dist kept
            %konstant --> get ver_dist with new slope
            %(tand(wiggled_slope)*hor_dist)
            src_ePot = g*tand(wiggled_slope).*hor_dist(j,i,:).*(slope<=0); %potential energy of source cell: (slope<=0)--> wiggle only downward slopes (with outflow)
            src_ePot = src_ePot + ePot(j,i,:).*(slope>=0); %include again upward slope --> taken from ePot which was not wiggled
            e_Kin = energy(j,i)+src_ePot-eFric(j,i,:); 
        else %after source cell --> remove potential energy
            e_Kin = energy(j,i)+ePot(j,i,:)*(1-en_decay)-eFric(j,i,:);
        end
        e_Kin = e_Kin(:)';
        %set v to v_max if larger than v_max; and calculate energy
        %again
        v = sqrt(2*e_Kin.*(e_Kin>0));
        v(v>v_max) = v_max;
        e_Kin = 0.5*(v).^2;
        %stop spreading if no energy available
        wgt_suscept(e_Kin<=0) = 0;

        %normalize wgt_susept again
        wgt_suscept = wgt_suscept./sum(wgt_suscept);

        %spread temporarly to check if values are smaller than delta_i
        %if yes --> set zero and redistribute by normalize again
        temp_spread = spread(j,i).*wgt_suscept;
        wgt_suscept(temp_spread<=delta_i) = 0;

        %normalize wgt_sucept again
        wgt_suscept = wgt_suscept./sum(wgt_suscept);

        %final spreading
        temp_spread = spread(j,i).*wgt_suscept;

        spread(j,i) = 0;

        %iteration through neighbours --> spread intensity
        %if it flows in an already active cell (also temporary active)
        %it sums up the intensity such that it is spread in the next
        %iteration and not lost. It takes the maximum of energy and
        %also the corresponding direction
        for c=1:8
            if temp_spread(c) > 0;
                %coordinates of neighbour

                J = j+shift_matrix(c,1); 
                I = i+shift_matrix(c,2);

                %spread_old = spread(j+shift_matrix(c,1),i+shift_matrix(c,2));
                %if temp_spread(c) > spread_old
                if temp_active_cells(J,I) || active_cells(J,I)
                    spread(J,I) = spread(J,I)+temp_spread(c);
                    %if new energy greater --> take new energy and its
                    %direction
                    old_energy = energy(J,I);
                    if e_Kin(c) > old_energy
                        energy(J,I) = e_Kin(c);
                        %energy(J,I) = e_Kin(c)*(1-en_decay);
                        direction(J,I) = c;
                    end
                else 
                    spread(J,I) = temp_spread(c);
                    energy(J,I) = e_Kin(c);
                    %energy(J,I) = e_Kin(c)*(1-en_decay);
                    direction(J,I) = c;
                end
                temp_active_cells(J,I) = 1;
                rel_prob(J,I) = rel_prob(J,I)+temp_spread(c);
                %save distance from source --> by adding distanc from
                %previous cell to inflow-cell (
                if cal_int
                    intensity(J,I) = intensity(j,i)* wgt_suscept(c)+cell_area(J,I);
                end
            end %end if spread > 0
        end %end interation through neighbours
        %active cell is treated --> set to zero
        active_cells(j,i)=0;
        direction(j,i) = 0;
        energy(j,i) = 0;
    %end %end if active cell
end %end iteration through active cells
active_cells(temp_active_cells == 1) = 1;
%%%%%for animation can be removed%%%%%
% if mod(k,3)==0
%     prob_plot = rel_prob;
%     prob_plot(prob_plot==0) = nan;
%     s.CData = log(prob_plot);
%     pause(1)
%     picnumber = picnumber+1;
%     fullpath = [path filesep 'pic' num2str(picnumber) '.tif'];
%     K = GRIDobj(lon,lat,log(prob_plot));
%     GRIDobj2geotiff(K,fullpath);
% end
% k
% k =k+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end %end when there are no active cells left


end