function climada_ls_test_flowParameters()

% Read in generated landslide sets (with intensity coordinates (lat/lon) 
% and elevation. With this information try out flowpath with different 
% parameters. Different analyses for master thesis. Not for the enduser.
% MODULE:
%   
% NAME:
%   
% PURPOSE:
%   
% CALLING SEQUENCE:
%   
% EXAMPLE:
%   
% INPUTS:
%   
% OPTIONAL INPUT PARAMETERS:
%  
% OUTPUTS:
%  
% MODIFICATION HISTORY:
% Thomas Rölli, thomasroelli@gmail.com, 20180306, init 
% Thomas Rölli, thomasroelli@gmail.com, 20180308, with animation --> gifs


%load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_hazard.mat')
load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_srtm1_hazard.mat')

%load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_centroids.mat')
load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_srtm1_centroids.mat')

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

%select subset of DEM
sublon = [60,250];
sublat = [110,300];
lon = lon(sublat(1):sublat(2),sublon(1):sublon(2));
lat = lat(sublat(1):sublat(2),sublon(1):sublon(2));
elevation = elevation(sublat(1):sublat(2),sublon(1):sublon(2));
intensity = intensity(sublat(1):sublat(2),sublon(1):sublon(2));
n_lon = numel(lon(1,:));
n_lat = numel(lat(:,1));

deg_km = 111.32;
dlat = abs(min(diff(lat(:,1)))); 
dlon = abs(min(diff(lon(1,:))));
dy = dlat*(deg_km * 1000);
dx = dlon*cosd(mean(lat(:,1)))*(deg_km * 1000); 

source = intensity(:,:,1);

elevation = deminpaint(elevation);
elevation = fillsinks(elevation);

[~,horDist,verDist] = climada_centroids_gradients(lon,lat,elevation);


%view([0,45])

%default parameters (exponent multipleflow = 25; max slide velocity = 8
%friction parameter phi = 18
% example: spreaded = climada_ls_flowpath(lon,lat,elevation,intensity,exponent,v_max,phi);
% exponent = 25;
% v_max = 8;
% phi = 18;


%plot 3*3*3 subplots with different exp,v-max,phi values
% for exp = exponent-1:exponent+1
%     for v = v_max-1:v_max+1
%         figure('units','normalized','outerposition',[0 0 1 1])
%         i=1;
%         for p = phi-1:phi+1
%             spreaded = climada_ls_flowpath(lon,lat,elevation,intensity,exp,v,p);
%             %figure('units','normalized','outerposition',[0 0 1 1])
%             subplot(1,3,i)
%             surf(lon,lat,elevation,double(spreaded(:,:,1)>0));
%             view([45,68])
%             title(['exp: ' num2str(exp) ', v_m: ' num2str(v) ', phi: ' num2str(p)]);
%             xlabel('lon');
%             ylabel('lat');
%             i=i+1;
%         end
%     end
% end

% k=9;
% figure('units','normalized','outerposition',[0 0 1 1])
% i=1;
% for exp = exponent-floor(k/2):exponent+floor(k/2)
%     spreaded = climada_ls_flowpath(lon,lat,elevation,intensity,exp,v_max,phi);
%     subplot(3,3,i)
%     surf(lon,lat,elevation,double(spreaded(:,:,1)>0));
%     view([45,68])
%     title(['exp: ' num2str(exp) ', v_m: ' num2str(v_max) ', phi: ' num2str(phi)]);
%     xlabel('lon');
%     ylabel('lat');
%     i=i+1;
% end
  

%select source cells(buffer region)
sel_source = zeros(size(elevation));
mask = zeros(size(elevation));

buf_m = 900; %bufferregion in meters in which no other slides are choosen--> prevents slides from flow over each other
imask = ceil(buf_m/dy);
jmask = ceil(buf_m/dx);
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

%%%%%%%%%%
%%%%GIF%%%
exp = 10;
v_max = 15;
phi = 50:-1:1;
pov = [45,68];

filename = 'anim_phi.gif';
%set model surface plot
h = figure('units','normalized','outerposition',[0 0 1 1]);
surf(lon,lat,elevation)
xlabel('longitude [\circ]','FontSize',20);
ylabel('latitude [\circ]','FontSize',20);
zlabel('elevation [m]','FontSize',20);
view(pov)
%zoom(1.3)
zl = zlim;
axis tight
zlim manual
set(gca,'nextplot','replacechildren')
firstpic = 1;
for e = 1:length(exp)
    for v = 1:length(v_max)
        for p = 1:length(phi)
            mult_flow = climada_ls_multipleflow(lon,lat,elevation,exp(e));
            spreaded = climada_ls_propagation(sel_source,mult_flow,horDist,verDist,v_max(v),phi(p));
            spreaded(spreaded>0)=1;
            spreaded(sel_source>0)=2;
            surf(lon,lat,elevation,double(spreaded(:,:,1)),'LineStyle','-') 
            view(pov)
            %zoom(1.3)
            title(['exp: ' num2str(exp(e)) ', v_{max}: ' num2str(v_max(v)) ', phi: ' num2str(phi(p))],'FontSize',20);
            drawnow 
            % Capture the plot as an image 
            frame = getframe(h); 
            im = frame2im(frame); 
            [imind,cm] = rgb2ind(im,256); 
            % Write to the GIF File 
            if firstpic == 1 
                imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
                firstpic = firstpic+1;
            else 
                imwrite(imind,cm,filename,'gif','DelayTime',0.1,'WriteMode','append'); 
            end 
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%compare old spread version and new%%%%

% exponent = 5;
% dH = 1;
% flat_areas = 1;
% v_max = 12;
% phi = 11;
% friction = 1;
% delta_i = 0.0001;
% %perWt = [1 1 1 1 1 1 1 1];
% perWt = [1 0.8 0.4 0 0 0 0.4 0.8];
% 
% mult_flow = climada_ls_multipleflow(lon,lat,elevation,exponent,dH,flat_areas);
% [~,hor_dist,ver_dist] = climada_centroids_gradients(lon,lat,elevation,dH);


%spread_old = climada_ls_spread(intensity(:,:,1),mult_flow,hor_dist,ver_dist,v_max,phi,friction);
% figure
% surface(spread_old)
% figure
% surface(log(spread_old))
% figure
% surface(spread_old>0)

%source_area = intensity(:,:,1);
% source_area = zeros(size(elevation));
% source_area(253,51) = 1;
% 
% spread_new = climada_ls_propagation(source_area,mult_flow,hor_dist,ver_dist,v_max,phi,delta_i,perWt);
% figure
% surface(elevation,spread_new)
% figure
% surface(elevation,log(spread_new))
% figure
% surface(elevation,spread_new>0)

%%%%%test new spread version spread_v2 with iteration through all events an

disp('hier')



end
