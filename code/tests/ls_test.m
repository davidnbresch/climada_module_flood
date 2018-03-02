function ls_test()
%function to test stuff... can be deleted afterwards

load('C:\Users\Simon Rölli\Desktop\data\z.mat','z');
load('C:\Users\Simon Rölli\Desktop\data\lon.mat','lon');
load('C:\Users\Simon Rölli\Desktop\data\lat.mat','lat');

exponent = 1.1;

[gradients,~,~] = climada_centroids_gradients(lon,lat,z);

outflow_gradients = gradients*-1;
outflow_gradients(outflow_gradients < 0) = 0; 

%%% calculate sum of all outflow cells
outflow_gradients_sum = sum(outflow_gradients.^exponent,3);
outflow_gradients_sum(outflow_gradients_sum==0) = 1; %prevent division by 0

%%% calculate multidirectional outflow proportion
% (tan(beta_i)^x/sum(tan(beta_i)^x(i= 1 to 8))
outflow_proportion = (outflow_gradients.^exponent)./outflow_gradients_sum;
inflow_proportion = circshift(outflow_proportion,4,3);
inflow_gradients = circshift(gradients,4,3);


total_field = outflow_gradients_sum*0;
field = outflow_gradients_sum*0+1;
inflow_temp = outflow_gradients*0;
shift_matrix = [-1 0;-1 -1;0 -1;1 -1;1 0;1 1;0 1;-1 1];

figure('units','normalized','outerposition',[0 0 1 1])
s = surf(lon,lat,z,total_field);
view([30,30,50])
colorbar
count = 0;
while sum(field(:))>0
    for i=1:8
        %inflow_temp(:,:,i) = circshift(field,shift_matrix(i,:)).*circshift(inflow_proportion(:,:,i),shift_matrix(i,:)).*(gradients(:,:,i)>0);
        inflow_temp(:,:,i) = circshift(field,shift_matrix(i,:)).*inflow_proportion(:,:,i).*(gradients(:,:,i)>0);
    end
    field = sum(inflow_temp,3);
    total_field = total_field + field;
    s.CData = total_field;
    pause(0.05)
    count = count+1
end


disp('hier')

end