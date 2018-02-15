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
%load('C:\Users\Simon Rölli\Desktop\climada\climada_data\centroids\_LS_Sarnen_centroids.mat')
load('C:\Users\Simon Rölli\Desktop\climada\climada_data\hazards\_LS_Sarnen_hazard.mat')

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
   [centroids,hazard] = climada_ls_testDEM;
else
   load('C:\Users\Simon Rölli\Desktop\climada\climada_data\centroids\_LS_Sarnen_centroids.mat') 
end

%get dimension of grid field from lon/lat coordinates
%and reshap needed vectors --> easier to handel in grid format than in
%vector; only possible for regular placed gridpoints
%also flip matrix up down to have high latitude in beginning of array
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

% mutiplying by 1 --> outflow should be positive; neighbours with inflow are set to zero
gradients = gradients*-1;
gradients(gradients < 0) = 0; 


%%% calculate sum of all outflow cells
gradients_sum = sum(gradients.^exponent,3);
gradients_sum(gradients_sum==0) = 1; %prevent division by 0

%%% calculate multidirectional outflow proportion
% (tan(beta_i)^x/sum(tan(beta_i)^x(i= 1 to 8))
mult_flow = (gradients.^exponent)./gradients_sum;


%shif matrix such that intensity is spread from center to neighbour-cells
%starting at 12 o'clock and proceeding clockwise
shift_matrix = [1 0;1 -1;0 -1;-1 -1;-1 0;-1 1;0 1;1 1]*-1;
max_runs=10;
g = 9.81; %acceleration of gravity
phi = 11; %empirical minimum travel angle, used for friction-calculation

%assessing flow path; starting from cells which are equal 1 in
%hazard.intensity
spread = double(zeros(n_lat,n_lon,hazard.event_count));
active_cells = logical(zeros(n_lat,n_lon)); %new for each event, save cells for next iteration
energy = double(zeros(n_lat,n_lon)); %new for each event, save energy of cells while propagating
eFric = g*horDist*tand(phi); %loss of energy by friction, for whole field in each direction
ePot = g*verDist*(-1); %gain of potential enenergy, for whole field in each direction

%%%%%%%without friction%%%%%%

for n_event=1:hazard.event_count %iteration through events
spread(:,:,n_event) = intensity(:,:,n_event);
active_cells = intensity(:,:,n_event);
while sum(sum(active_cells))>0 %iteration through number of runs --> implement friction here 
temp_active_cells = logical(zeros(n_lat,n_lon));
for j=1:n_lat %iteration through rows
    for i=1:n_lon %iteration through colums
        if active_cells(j,i)
            for c=1:8 %iteration through shift matrix
                if mult_flow(j,i,c) > 0 
                    %spread value of center to neighbour cell only if
                    %new value greater than old
                    if spread(j,i,n_event)*mult_flow(j,i,c) > spread(j+shift_matrix(c,1),i+shift_matrix(c,2),n_event)
                        %intensity is spread according to its outflow
                        %propotion
                        spread(j+shift_matrix(c,1),i+shift_matrix(c,2),n_event)=spread(j,i,n_event)*mult_flow(j,i,c);
                        temp_active_cells(j+shift_matrix(c,1),i+shift_matrix(c,2))=1;
                    end
                    %test = zeros(n_lon,n_lat);
                    %test(spread(:,:,1)>0) =1;
                    %test = double(active_cells);
                    %surf(lon,lat,elevation,test);
                    %view(2);
                    %disp([j i]);
                    %disp('');
                end
            end %end interation through shift matrix
            active_cells(j,i)=0;
        end
    end %end interation through columns
end %end interation through rows
active_cells(temp_active_cells == 1) = 1;
end %end interation through outflow distance
end %end interation through events


%%%%%%%with friction%%%%%%
spreadFri = double(zeros(n_lat,n_lon,hazard.event_count));

for n_event=1:hazard.event_count %iteration through events
spreadFri(:,:,n_event) = intensity(:,:,n_event);
active_cells = intensity(:,:,n_event);
energy(energy~=0) = 0; %set energy to zero for new event
while sum(sum(active_cells))>0 %iteration through number of runs --> implement friction here 
temp_active_cells = logical(zeros(n_lat,n_lon));    
for j=1:n_lat %iteration through rows
    for i=1:n_lon %iteration through colums
        if active_cells(j,i)
            for c=1:8 %iteration through shift matrix
                if mult_flow(j,i,c) > 0
                    %kin. energy to corresponding neighbour
                    eKin = energy(j,i)+ePot(j,i,c)-eFric(j,i,c);
                    spread_old = spreadFri(j+shift_matrix(c,1),i+shift_matrix(c,2),n_event);
                    %calculation of new intensity when spreading
                    spread_new = spreadFri(j,i,n_event)*mult_flow(j,i,c);
                    %spread value of center to neighbour cell only if
                    %new value greater than old and there is kinetic energy
                    if (spread_new > spread_old) & (eKin > 0)
                        spreadFri(j+shift_matrix(c,1),i+shift_matrix(c,2),n_event)=spread_new;
                        energy(j+shift_matrix(c,1),i+shift_matrix(c,2)) = eKin;
                        temp_active_cells(j+shift_matrix(c,1),i+shift_matrix(c,2))=1;
                    end
                end
            end %end interation through shift matrix
            active_cells(j,i)=0;
        end
    end %end interation through columns
end %end interation through rows
active_cells(temp_active_cells == 1) = 1;
end %end interation through outflow distance
end %end interation through events

figure
surf(lon,lat,elevation,spread(:,:,1));
figure
surf(lon,lat,elevation,spreadFri(:,:,1));

test = zeros(n_lon,n_lat);
test(spreadFri(:,:,1)>0) =1;
figure
surf(lon,lat,elevation,test)
disp('hier');












end
