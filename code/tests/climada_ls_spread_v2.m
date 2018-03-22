function spread = climada_ls_spread_v2(source_area,mult_flow,... 
    hor_dist,ver_dist,v_max,phi,delta_i,perWt)

% Computes the flow path according to the multiple flow algorithm
% (according to Holmgren 1994). The flow distance is taken into account by
% a simplified friction model (see Horton et al. 2013).
% MODULE:
%   flood
% NAME:
%   climada_ls_spread
% PURPOSE:
%   
% CALLING SEQUENCE:
%   climada_ls_hazard_sets
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
%   friction: 1/0 include friction/no friction while spreading
%   v_max: describes maximal possible velocity of slide. If velocity is 
%          exceeded v_max is taken. Should keep energy amounts within reasonal values
%          and therefore prevent improbalbe runout distances
%   phi: angle of reach, angle of the line connnecting the source area to
%        the most distant point reached by the slide, along its path.
%        Factor controlls maximum possible runout distance
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
% OPTIONAL INPUT PARAMETERS:
% 
% OUTPUTS:
%   spread: matrix (lon/lat) where the intensity of the landslide is spread
%   according to the multiple flow algorithm.
%  
% MODIFICATION HISTORY:
% Thomas R�lli, thomasroelli@gmail.com, 20180219, init
% Thomas R�lli, thomasroelli@gmail.com, 20180227, changed shift_matrix
% Thomas R�lli, thomasroelli@gmail.com, 20180305, implement v_max of slide
% Thomas R�lli, thomasroelli@gmail.com, 20180317, init v2
% Thomas R�lli, thomasroelli@gmail.com, 20180319, changed structure, add
%  persistence function, add intensity threshold (delta_i)


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

% PARAMETERS 
if isempty(source_area); return; end
if isempty(mult_flow); return; end
if isempty(hor_dist); return; end
if isempty(ver_dist); return; end
if isempty(v_max); v_max = 8; end
if isempty(phi); phi = 18; end %empirical minimum travel angle, used for friction-calculation
if isempty(delta_i); delta_i = 0.0001; end
if isempty(perWt); perWt = [1 0.8 0.4 0 0 0 0.4 0.8]; end

if numel(perWt) ~= 8; return; end

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


%direction-matrix --> saves direction of flow (1-8)in each iteration
direction = source_area*0; %will save flow-direction in each iteration

while sum(sum(active_cells))>0 %iteration through number of runs  
temp_active_cells = temp_active_cells*0;
for j=1:n_lat %iteration through rows
    for i=1:n_lon %iteration through colums
        if active_cells(j,i)
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
            e_Kin = energy(j,i)+ePot(j,i,:)-eFric(j,i,:);
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
            
            %write in spread-matrix if value is larger than old value
            for c=1:8
                if temp_spread(c) > 0;
                    spread_old = spread(j+shift_matrix(c,1),i+shift_matrix(c,2));
                    if temp_spread(c) > spread_old
                        spread(j+shift_matrix(c,1),i+shift_matrix(c,2)) = temp_spread(c);
                        energy(j+shift_matrix(c,1),i+shift_matrix(c,2)) = e_Kin(c);
                        direction(j+shift_matrix(c,1),i+shift_matrix(c,2)) = c;
                        temp_active_cells(j+shift_matrix(c,1),i+shift_matrix(c,2)) = 1;
                    end
                end
            end %end interation through neighbours
            %active cell is treated --> set to zero
            active_cells(j,i)=0;
            direction(j,i) = 0;
        end %end if active cell
    end %end interation through columns
end %end interation through rows
active_cells(temp_active_cells == 1) = 1;
end %end when there are no active cells left


disp('hier')


end