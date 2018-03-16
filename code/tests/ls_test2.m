function ls_test2()
%function to test stuff... can be deleted afterwards

%load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_centroids.mat')
load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_srtm1_centroids.mat')

%load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_hazard.mat')
%load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_srtm1_hazard.mat')

n_lon = numel(unique(centroids.lon));
n_lat = numel(unique(centroids.lat));
lon = reshape(centroids.lon,n_lat,n_lon);
lat = reshape(centroids.lat,n_lat,n_lon);
elevation = reshape(centroids.elevation_m,n_lat,n_lon);

dH = 0;
flat = 1;
g = 9.81;
phi = 12;
v_max = 8;
friction = 1;
exponent = 1;
delta_i = 0.0001;

mult_flow = climada_ls_multipleflow(lon,lat,elevation,exponent,dH,flat);

[~,hor_dist,ver_dist] = climada_centroids_gradients(lon,lat,elevation,dH);

index = [200 50];
source_area = zeros(n_lat,n_lon);
source_area(index(1),index(2)) = 1;

%shif matrix such that intensity is spread from center to neighbour-cells
%starting at 12 o'clock and proceeding clockwise
%shift_matrix = [1 0;1 -1;0 -1;-1 -1;-1 0;-1 1;0 1;1 1]*-1;
shift_matrix = [-1 0;-1 -1;0 -1;1 -1;1 0;1 1;0 1;-1 1]*-1;

%assessing flow path; starting from cells which are equal 1 in
spread = source_area;
active_cells = source_area; %save active cells for next iteration
temp_active_cells = active_cells*0;
energy = double(zeros(n_lat,n_lon)); %save energy of cells while propagating
eFric = g*hor_dist*tand(phi); %loss of energy by friction, for whole field in each direction
ePot = g*ver_dist*(-1); %gain of potential enenergy, for whole field in each direction

%persistence weighting matrix and direction of flow; first number is when there is no
%direction --> for source areas
perWt = [1 0.8 0.4 0 0 0 0.4 0.8];
perSum = sum(perWt);
direction = source_area*0;

spreaded = spread>0;

figure('units','normalized','outerposition',[0 0 1 1])
fov = 49;
s = surface(elevation(index(1)-fov:index(1)+fov,index(2)-fov:index(2)+fov),spreaded(index(1)-fov:index(1)+fov,index(2)-fov:index(2)+fov));
colorbar
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
%weitermachen: susceptibility wird berechnet aus produkt mit persistence
%und mult_flow --> zuerst per*multflow (mit circshift) --> normieren mit
%summe nicht vergessen.
%dann susceptiblity verteilen gemäss overallausflussanteil (per und mult
%flow inbegriffen)
%danach schauen ob zellen ereicht werden (energy delta_i) wenn nicht -->
%eintrag löschen und neu overall normieren dass 1 ergibt.
%dann in schleife durch 8 richtungen (wenn overall susceptiblity > 0) --> neue susceptibility berechnen (nur
%wenn alt kleinere suceptiblity als neu; alle werte einfügen und
%temp_active cell activieren
%am ende active_cell deaktiviren und tem_active_cells einfügen




%spreading for single point
while sum(sum(active_cells))>0 %iteration through number of runs  
temp_active_cells = temp_active_cells*0;
for j=1:n_lat %iteration through rows
    for i=1:n_lon %iteration through colums
        if active_cells(j,i)
            suceptiblity = zeros(8,1);
            temp_mult_flow = mult_flow(j,i,:);
            temp_mult_flow = temp_mult_flow(:);
            if direction(j,i) ~= 0 
                shifted_perWt = circshift(perWt,direction(j,i)-1);
            else
                shifted_perWt = ones(8,1);
            end
                
            shifted_perWt = circshift(perWt,direction(j,i))
            %temp_mult_flow = mult_flow(j,i,:);
            
            %check if next cells are reached according to threshold delta_i
            %and/or kinetic energy; if not intensity is redistributed to
            %other neighbours which are reached
%             temp_mult_flow = mult_flow(j,i,:);
%             spread_new = spread(j,i)*temp_mult_flow;
%             gt_delta_i = spread_new > delta_i;
%             
%             %check energy to its neighbour
%             e_Kin = energy(j,i)+ePot(j,i,:)-eFric(j,i,:);
%             v = sqrt(2*e_Kin.*(e_Kin>0));
%             %if velocity to high --> set to v_max
%             v(v>v_max) = v_max;
%             e_Kin = 0.5*(v).^2; 
            
            
            
            
            for c=1:8 %iteration through shift matrix
                if mult_flow(j,i,c) > 0
                    %calculation of new intensity when spreading
                    spread_new = spread(j,i)*mult_flow(j,i,c);
                    spread_old = spread(j+shift_matrix(c,1),i+shift_matrix(c,2));
                    %spread value of center to neighbour cell only if
                    %new value greater than old
                    if spread_new > spread_old
                    switch friction
                        case true
                            %kin. energy if spread to corresponding neighbour
                            eKin = energy(j,i)+ePot(j,i,c)-eFric(j,i,c);
                            %spread only when there is kinetic energy
                            if eKin > 0
                                %calculate vi; if vi > vmax --> set vi to vmax
                                v_i = sqrt(2*eKin);
                                v_i = min(v_i,v_max);
                                eKin = 0.5*(v_i)^2;
                                spread(j+shift_matrix(c,1),i+shift_matrix(c,2))=spread_new;
                                energy(j+shift_matrix(c,1),i+shift_matrix(c,2)) = eKin;
                                direction(j+shift_matrix(c,1),i+shift_matrix(c,2))=c;
                                temp_active_cells(j+shift_matrix(c,1),i+shift_matrix(c,2))=1; 
                            end
                        case false
                            spread(j+shift_matrix(c,1),i+shift_matrix(c,2))=spread_new;
                            direction(j+shift_matrix(c,1),i+shift_matrix(c,2))=c;
                            temp_active_cells(j+shift_matrix(c,1),i+shift_matrix(c,2))=1;
                    end
                    end
                end
            end %end interation through shift matrix
            active_cells(j,i)=0;
            direction(j,i) = 0;
        end
    end %end interation through columns
end %end interation through rows
active_cells(temp_active_cells == 1) = 1;
active_cells(j,i)=0;
spreaded = double(spread>0);
s.CData = spreaded(index(1)-fov:index(1)+fov,index(2)-fov:index(2)+fov);
pause(0.5)
end %end when there are no active cells left




disp('hier')


end