function ls_test()
%function to test stuff... can be deleted afterwards

% load('C:\Users\Simon Rölli\Desktop\data\z.mat','z');
% load('C:\Users\Simon Rölli\Desktop\data\lon.mat','lon');
% load('C:\Users\Simon Rölli\Desktop\data\lat.mat','lat');

load('C:\Users\Simon Rölli\Desktop\data\centroids_hazards_nospread\_LS_Sarnen_srtm1_centroids.mat')

%get gridded datasets
n_lon = numel(unique(centroids.lon));
n_lat = numel(unique(centroids.lat));
lon = reshape(centroids.lon,n_lat,n_lon);
lat = reshape(centroids.lat,n_lat,n_lon);
z = reshape(centroids.elevation_m,n_lat,n_lon);

z = deminpaint(z);

%z = fillsinks(z);

mult_flow = climada_ls_multipleflow(lon,lat,z,1.1,0);


%   DEM = GRIDobj(lon,lat,flipud(z));
% % 
% %  DEM = fillsinks(DEM);
% %  
% [I,SILLS] = identifyflats(DEM);
%  
%  fl = FLOWobj(DEM,'multi');
% 
% %fl = FLOWobj(DEM,'multi');
% 
% DEM = imposemin(fl,DEM);
% 
% z = DEM.Z;



% exponent = 1.1;
% 
% [gradients,~,~] = climada_centroids_gradients(lon,lat,z);
% 
% outflow_gradients = gradients*-1;
% outflow_gradients(outflow_gradients < 0) = 0; 
% 
% %%% calculate sum of all outflow cells
% outflow_gradients_sum = sum(outflow_gradients.^exponent,3);
% outflow_gradients_sum(outflow_gradients_sum==0) = 1; %prevent division by 0
% 
% %%% calculate multidirectional outflow proportion
% % (tan(beta_i)^x/sum(tan(beta_i)^x(i= 1 to 8))
% outflow_proportion = (outflow_gradients.^exponent)./outflow_gradients_sum;
inflow_proportion = circshift(mult_flow,4,3);



total_field = z*0;
field = z*0+1;
inflow_temp = z*0;
shift_matrix = [-1 0;-1 -1;0 -1;1 -1;1 0;1 1;0 1;-1 1];

%figure('units','normalized','outerposition',[0 0 1 1])
%s = surf(lon,lat,z,total_field);
%view([30,30,50])
%colorbar
count = 0;
while sum(field(:),'omitnan')>0
    for i=1:8
        inflow_temp(:,:,i) = circshift(field,shift_matrix(i,:)).*circshift(inflow_proportion(:,:,i),shift_matrix(i,:));
        %inflow_temp(:,:,i) = (field,shift_matrix(i,:)).inflow_proportion(:,:,i);
    end
    field = sum(inflow_temp,3);
    total_field = total_field + field;
    %s.CData = total_field;
    %pause(0.05)
    %count = count+1
end

figure
surface(z,log(total_field),'LineStyle','none')
disp('hier')

end