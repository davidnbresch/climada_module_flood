function ls_test2()
%function to test stuff... can be deleted afterwards

load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_hazard.mat')
%load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_srtm1_hazard.mat')

load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_centroids.mat')
%load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_srtm1_centroids.mat')

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

deg_km = 111.32;
dlat = abs(min(diff(lat(:,1)))); 
dlon = abs(min(diff(lon(1,:))));
dy = dlat*(deg_km * 1000);
dx = dlon*cosd(mean(lat(:,1)))*(deg_km * 1000);

source = intensity(:,:,1);

elevation = deminpaint(elevation);
elevation = fillsinks(elevation);

exp = 10;
v_max = 12;
phi = 10;

[~,horDist,verDist] = climada_centroids_gradients(lon,lat,elevation);
mult_flow = climada_ls_multipleflow(lon,lat,elevation,exp);

cell_area = climada_centroids_area(lon,lat,elevation);

source = source*0;
source(80,15)=1;
[spreaded,intensity] = climada_ls_propagation(source,mult_flow,horDist,verDist,v_max,phi,'','',1,cell_area);

figure
surface(lon,lat,intensity)

%select source cells(buffer region)
sel_source = zeros(size(elevation));
mask = zeros(size(elevation));
tot_spreaded = zeros(size(elevation));



buf_m = 1000; %bufferregion in meters in which no other slides are choosen--> prevents slides from flow over each other
imask = ceil(buf_m/dy);
jmask = ceil(buf_m/dx);
tic
while sum(source(:)) > 0
    for i=1:n_lat
        for j=1:n_lon
            if (source(i,j) == 1) && (mask(i,j) ~= 1)
                sel_source(i,j) = 1;

                %set values in mask within range to 1 --> no slides in this
                %area

                Imin=i-imask;Imax=i+imask;
                Jmin=j-jmask;Jmax=j+jmask;
                if (Imin < 1) Imin=1; end
                if (Imax > n_lat) Imax=n_lat; end
                if (Jmin < 1) Jmin=1; end
                if (Jmax > n_lon) Jmax=n_lon; end
                mask(Imin:Imax,Jmin:Jmax) = 1;
                %remove selected cell
                source(i,j) = 0;
            end 
        end %interation through colums
    end %iteration through rows
    spreaded = climada_ls_propagation(sel_source,mult_flow,horDist,verDist,v_max,phi);
    tot_spreaded = max(tot_spreaded,spreaded);

    mask(:) = 0;
    sel_source(:) = 0;
end
toc

tic
tot_spreaded_single = zeros(size(elevation));
source = intensity(:,:,1);
for i=1:n_lat
    for j=1:n_lon
        if source(i,j) == 1
            scr = zeros(size(elevation));
            scr(i,j) = 1;
            spreaded = climada_ls_propagation(scr,mult_flow,horDist,verDist,v_max,phi);
            %tot_spreaded_single = max(tot_spreaded_single,spreaded);
        end
    end
end
toc          

%test of old (slow version --> climada_ls_propagation_slow has been removed
%--> see gitHub older version of climada_ls_propagation for old/slow
%version
% sel_source = intensity(:,:,1);
% sum(sum(sel_source))
% k=0;
% tic
% tot_spreaded_slow = zeros(size(elevation));
% for i=1:n_lat
%     for j=1:n_lon
%         if sel_source(i,j) == 1
%             k=k+1;
%             if mod(k,100) == 0
%                k
%             end
%             scr = zeros(size(elevation));
%             scr(i,j) = 1;
%             spreaded = climada_ls_propagation_slow(scr,mult_flow,horDist,verDist,v_max,phi);
%             tot_spreaded_slow = max(tot_spreaded_slow,spreaded);
%         end
%     end
% end
% toc


end