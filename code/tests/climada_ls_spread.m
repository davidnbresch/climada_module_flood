function [spread,spreadFri] = climada_ls_spread(source_area,mult_flow,... 
    hor_dist,ver_dist,friction)

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
% OPTIONAL INPUT PARAMETERS:
% 
% OUTPUTS:
%   spread: matrix (lon/lat) where the intensity of the landslide is spread
%   according to the multiple flow algorithm.
%  
% MODIFICATION HISTORY:
% Thomas Rölli, thomasroelli@gmail.com, 20180219, init
% Thomas Rölli, thomasroelli@gmail.com, 20180227, changed shift_matrix



global climada_global
if ~climada_init_vars, return; end

% check arguments
if ~exist('source_area', 'var'), source_area = []; end
if ~exist('mult_flow', 'var'), mult_flow = []; end
if ~exist('hor_dist', 'var'), hor_dist = []; end
if ~exist('ver_dist', 'var'), ver_dist = []; end
if ~exist('friction', 'var'), friction = []; end

% PARAMETERS 
if isempty(source_area); return; end
if isempty(mult_flow); return; end
if isempty(hor_dist); return; end
if isempty(ver_dist); return; end
if isempty(friction); friction = true; end
%for calculations
g = 9.81; %acceleration of gravity
phi = 11; %empirical minimum travel angle, used for friction-calculation

%shif matrix such that intensity is spread from center to neighbour-cells
%starting at 12 o'clock and proceeding clockwise
%shift_matrix = [1 0;1 -1;0 -1;-1 -1;-1 0;-1 1;0 1;1 1]*-1;
shift_matrix = [-1 0;-1 -1;0 -1;1 -1;1 0;1 1;0 1;-1 1]*-1;

%get dimension of study area
n_lat = numel(source_area(:,1,1));
n_lon = numel(source_area(1,:,1));
n_event = numel(source_area(1,1,:));

%assessing flow path; starting from cells which are equal 1 in
spread = double(zeros(n_lat,n_lon,n_event));
active_cells = logical(zeros(n_lat,n_lon)); %new for each event, save cells for next iteration
energy = double(zeros(n_lat,n_lon)); %new for each event, save energy of cells while propagating
eFric = g*hor_dist*tand(phi); %loss of energy by friction, for whole field in each direction
ePot = g*ver_dist*(-1); %gain of potential enenergy, for whole field in each direction



% %%%%%%%without friction%%%%%%
% 
% for event_count=1:n_event %iteration through events
% spread(:,:,event_count) = source_area(:,:,event_count);
% active_cells = source_area(:,:,event_count);
% while sum(sum(active_cells))>0 %iteration through number of runs --> implement friction here 
% temp_active_cells = logical(zeros(n_lat,n_lon));
% for j=1:n_lat %iteration through rows
%     for i=1:n_lon %iteration through colums
%         if active_cells(j,i)
%             for c=1:8 %iteration through shift matrix
%                 if mult_flow(j,i,c) > 0 
%                     %spread value of center to neighbour cell only if
%                     %new value greater than old
%                     if spread(j,i,event_count)*mult_flow(j,i,c) > spread(j+shift_matrix(c,1),i+shift_matrix(c,2),event_count)
%                         %intensity is spread according to its outflow
%                         %propotion
%                         spread(j+shift_matrix(c,1),i+shift_matrix(c,2),event_count)=spread(j,i,event_count)*mult_flow(j,i,c);
%                         temp_active_cells(j+shift_matrix(c,1),i+shift_matrix(c,2))=1;
%                     end
%                 end
%             end %end interation through shift matrix
%             active_cells(j,i)=0;
%         end
%     end %end interation through columns
% end %end interation through rows
% active_cells(temp_active_cells == 1) = 1;
% end %end interation through outflow distance
% end %end interation through events


%%%%%%%with friction%%%%%%
spread = double(zeros(n_lat,n_lon,n_event));

for event_count=1:n_event %iteration through events
spread(:,:,event_count) = source_area(:,:,event_count);
active_cells = source_area(:,:,event_count);
energy(energy~=0) = 0; %set energy to zero for new event
while sum(sum(active_cells))>0 %iteration through number of runs --> implement friction here 
temp_active_cells = logical(zeros(n_lat,n_lon));    
for j=1:n_lat %iteration through rows
    for i=1:n_lon %iteration through colums
        if active_cells(j,i)
            for c=1:8 %iteration through shift matrix
                if mult_flow(j,i,c) > 0
                    %calculation of new intensity when spreading
                    spread_new = spread(j,i,event_count)*mult_flow(j,i,c);
                    spread_old = spread(j+shift_matrix(c,1),i+shift_matrix(c,2),event_count);
                    %spread value of center to neighbour cell only if
                    %new value greater than old
                    if spread_new > spread_old
                    switch friction
                        case true
                            %kin. energy to corresponding neighbour
                            eKin = energy(j,i)+ePot(j,i,c)-eFric(j,i,c);
                            %spread only when there is kinetic energy
                            if eKin > 0
                                spread(j+shift_matrix(c,1),i+shift_matrix(c,2),event_count)=spread_new;
                                energy(j+shift_matrix(c,1),i+shift_matrix(c,2)) = eKin;
                                temp_active_cells(j+shift_matrix(c,1),i+shift_matrix(c,2))=1; 
                            end
                        case false
                            spread(j+shift_matrix(c,1),i+shift_matrix(c,2),event_count)=spread_new;
                            temp_active_cells(j+shift_matrix(c,1),i+shift_matrix(c,2))=1;
                    end
                    end
                end
            end %end interation through shift matrix
            active_cells(j,i)=0;
        end
    end %end interation through columns
end %end interation through rows
active_cells(temp_active_cells == 1) = 1;
end %end when there are no active cells left
end %end interation through events

end