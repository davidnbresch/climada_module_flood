function [f_plot,f_con] = climada_calibration_plotRMSE(res,res2,res3,res_,res2_,res3_)

% Ploting of the RMSE for all shapefiles defined in files
% 
% MODULE:
%   flood
% NAME:
%   climada_calibration_plotRMSE
% PURPOSE:
%   Plots the results of climada_calibration_calcRMSE in subplots and vmax
%   vs phi plot
% CALLING SEQUENCE:
%   climada_calibration_plotRMSE(res,rmv,res2,res3)
% EXAMPLE:
%   climada_ls_multipleflow(lon,lat,elevation,25,1)
% INPUTS: 
%   files:   structure (files = dir(['path' '/*.shp']) which contains filepaths 
%            of shapefiles (generated by scriptFPCalibration, containing model-output)
%            with different flow path parameters sets (need to be defined in
%            filename)
%            shapefile need to include fields: length, area
%   subS:    original shapefile (generated by climada_calibration_subS2raster)
%            whith observed length and area fields. Need to have same
%            number of slides as model-output shape files
%       
% OPTIONAL INPUT PARAMETERS:
%   rmv:     in
%   fig:     plot figures y/n 1/0
% OUTPUTS:
%   res:    
% MODIFICATION HISTORY:
% Thomas R�lli, thomasroelli@gmail.com, 20180523, init
% Thomas R�lli, thomasroelli@gmail.com, 20180604, add plots

%check arguments
if ~exist('res'), return; end
if ~exist('res2'), res2 = []; end
if ~exist('res3'), res3 = []; end
if ~exist('res_'), res_ = []; end
if ~exist('res2_'), res2_ = []; end
if ~exist('res3_'), res3_ = []; end

%read in data
vmax = [res.vmax];
vmax_unique = unique(vmax);
phi = [res.phi];
phi_unique = unique(phi);
rmse = [res.rmse_lgt];

if ~isempty(res2)
    vmax2 = [res2.vmax];
    phi2 = [res2.phi];
    rmse2 = [res2.rmse_lgt];
end

if ~isempty(res3)
    vmax3 = [res3.vmax];
    phi3 = [res3.phi];
    rmse3 = [res3.rmse_lgt];
end

if ~isempty(res_)
    vmax_ = [res_.vmax];
    phi_ = [res_.phi];
    rmse_ = [res_.rmse_lgt];
end

if ~isempty(res2_)
    vmax2_ = [res2_.vmax];
    phi2_ = [res2_.phi];
    rmse2_ = [res2_.rmse_lgt];
end

if ~isempty(res3_)
    vmax3_ = [res3_.vmax];
    phi3_ = [res3_.phi];
    rmse3_ = [res3_.rmse_lgt];
end

%%%%plot for phi

%figure with subplots
f_plot = figure('units','normalized','outerposition',[0 0 1 1]);
dim = [3 4];%dimension of subplot
subplot1(dim(1),dim(2),'Gap',[0.01 0.01],'XTickL','Margin','YTickL','Margin','FontS',10)

%create subplots
for i=1:numel(vmax_unique)
    %extract needed data for subplot
    idx_vmax = find(vmax == vmax_unique(i));
    phi_plot = phi(idx_vmax);
    rmse_plot = rmse(idx_vmax);
    %sort values according to phi
    [phi_plot,idx_sort] = sort(phi_plot);
    rmse_plot = rmse_plot(idx_sort);
    %plot subplot
    subplot1(i)
    p1=plot(phi_plot,rmse_plot,'-','color',[0, 0.4470, 0.7410],'LineWidth',2);
    [~,idx_min] = min(rmse_plot);
    grid on
    text(0.5,0.9,['VMAX = ' num2str(vmax_unique(i))],'FontSize',10,'FontWeight','bold',...
        'Units','normalized','HorizontalAlignment','center')
    %label of axis only at boarder
    if floor((i-1)/dim(2))+1==dim(1), xlabel('PHI [\circ]','FontWeight','bold'); end
    if mod(i-1,dim(2))==0, ylabel('RMSE','FontWeight','bold'); end
    xlim([min(phi_unique) max(phi_unique)])
    ylim([0 500])  
    if ~isempty(res2)
        phi_plot2 = phi2(idx_vmax);
        rmse_plot2 = rmse2(idx_vmax);
        [phi_plot2,idx_sort] = sort(phi_plot2);
        rmse_plot2 = rmse_plot2(idx_sort);
        p2=plot(phi_plot2,rmse_plot2,'-','color',[0.8500, 0.3250, 0.0980],'LineWidth',2);
        [~,idx_min2] = min(rmse_plot2);
        
    end
    if ~isempty(res3)
        phi_plot3 = phi3(idx_vmax);
        rmse_plot3 = rmse3(idx_vmax);
        [phi_plot3,idx_sort] = sort(phi_plot3);
        rmse_plot3 = rmse_plot3(idx_sort);
        p3=plot(phi_plot3,rmse_plot3,'-','color',[0.4660, 0.6740, 0.1880],'LineWidth',2);
        [~,idx_min3] = min(rmse_plot3);
    end
    if ~isempty(res_)
        phi_plot_ = phi_(idx_vmax);
        rmse_plot_ = rmse_(idx_vmax);
        [phi_plot_,idx_sort] = sort(phi_plot_);
        rmse_plot_ = rmse_plot_(idx_sort);
        plot(phi_plot_,rmse_plot_,'-.','color',[0, 0.4470, 0.7410],'LineWidth',1)
        [~,idx_min_] = min(rmse_plot_);
    end
    if ~isempty(res2_)
        phi_plot2_ = phi2_(idx_vmax);
        rmse_plot2_ = rmse2_(idx_vmax);
        [phi_plot2_,idx_sort] = sort(phi_plot2_);
        rmse_plot2_ = rmse_plot2_(idx_sort);
        plot(phi_plot2_,rmse_plot2_,'-.','color',[0.8500, 0.3250, 0.0980],'LineWidth',1)
        [~,idx_min2_] = min(rmse_plot2_);
    end
    if ~isempty(res3_)
        phi_plot3_ = phi3_(idx_vmax);
        rmse_plot3_ = rmse3_(idx_vmax);
        [phi_plot3_,idx_sort] = sort(phi_plot3_);
        rmse_plot3_ = rmse_plot3_(idx_sort);
        plot(phi_plot3_,rmse_plot3_,'-.','color',[0.4660, 0.6740, 0.1880],'LineWidth',1)
        [~,idx_min3_] = min(rmse_plot3_);
    end
    
    %plot minimum
    plot(phi_plot(idx_min),rmse_plot(idx_min),'o','MarkerSize',7,'color',[0, 0.4470, 0.7410],'MarkerFaceColor',[0, 0.4470, 0.7410])
    if ~isempty(res2)
        plot(phi_plot2(idx_min2),rmse_plot2(idx_min2),'s','MarkerSize',6,'color',[0.8500, 0.3250, 0.0980],'MarkerFaceColor',[0.8500, 0.3250, 0.0980])
    end
    if ~isempty(res3)
        plot(phi_plot3(idx_min3),rmse_plot3(idx_min3),'d','MarkerSize',5,'color',[0.4660, 0.6740, 0.1880],'MarkerFaceColor',[0.4660, 0.6740, 0.1880])
    end
    if ~isempty(res_)
        plot(phi_plot_(idx_min_),rmse_plot_(idx_min_),'o','MarkerSize',7,'color',[0, 0.4470, 0.7410],'MarkerFaceColor',[0, 0.4470, 0.7410])
    end
    if ~isempty(res2_)
        plot(phi_plot2_(idx_min2_),rmse_plot2_(idx_min2_),'s','MarkerSize',6,'color',[0.8500, 0.3250, 0.0980],'MarkerFaceColor',[0.8500, 0.3250, 0.0980])
    end
    if ~isempty(res3_)
        plot(phi_plot3_(idx_min3_),rmse_plot3_(idx_min3_),'d','MarkerSize',5,'color',[0.4660, 0.6740, 0.1880],'MarkerFaceColor',[0.4660, 0.6740, 0.1880])
    end
end

legend([p1,p2,p3],{'ALTI3D','SRTM1','SRTM3'})


%%%plot vmax vs phi
v = 1100:-100:100; %define contour classes
%define colormap
cmap = jet(numel(v)+1);
cmap = cmap(3:end,:);
xtic = min(phi):5:max(phi);
ytic = [1 5 10 12];

f_con = figure('units','normalized','outerposition',[0 0 1 1]);
subplot1(3,1,'Gap',[0.01 0.02],'XTickL','Margin','YTickL','Margin','FontS',10)

[VMAX,PHI] = meshgrid(linspace(min(phi),max(phi),200), linspace(min(vmax),max(vmax),200));
subplot1(1)
contourf(VMAX,PHI,griddata(phi,vmax,rmse,VMAX,PHI),v)
hold on
%plot grid
for i=1:numel(ytic)-2,plot([xtic(1),xtic(end)],[ytic(i+1),ytic(i+1)],'-.','Color','black');end
for i=1:numel(xtic)-2,plot([xtic(i+1),xtic(i+1)],[ytic(1),ytic(end)],'-.','Color','black');end
%plot frame
plot([xtic(1),xtic(end),xtic(end),xtic(1),xtic(1)],[ytic(1),ytic(1),ytic(end),ytic(end),ytic(1)],...
    '-','LineWidth',2,'color','black')
hold off
set(gca,'fontsize', 12)
colormap(cmap)
caxis([min(v) max(v)])
xticks(xtic)
yticks(ytic)
text(0.98,0.90,'ALTI3D','FontSize',15,'FontWeight','bold',...
        'Units','normalized','HorizontalAlignment','right')
axis equal
box on



if ~isempty(res2)
    [VMAX,PHI] = meshgrid(linspace(min(phi2),max(phi2),200), linspace(min(vmax2),max(vmax2),200));
    subplot1(2)
    contourf(VMAX,PHI,griddata(phi2,vmax2,rmse2,VMAX,PHI),v)
    %grid line in xdircetion and ydirection
    hold on
    for i=1:numel(ytic)-2,plot([xtic(1),xtic(end)],[ytic(i+1),ytic(i+1)],'-.','Color','black');end
    for i=1:numel(xtic)-2,plot([xtic(i+1),xtic(i+1)],[ytic(1),ytic(end)],'-.','Color','black');end
    plot([xtic(1),xtic(end),xtic(end),xtic(1),xtic(1)],[ytic(1),ytic(1),ytic(end),ytic(end),ytic(1)],...
    '-','LineWidth',2,'color','black')
    hold off
    set(gca,'fontsize', 12)
    ylabel('VMAX [m/s]','FontWeight','bold')
    caxis([min(v) max(v)])
    xticks(xtic)
    yticks(ytic)
    text(0.98,0.90,'SRTM1','FontSize',15,'FontWeight','bold',...
        'Units','normalized','HorizontalAlignment','right')
    axis equal
end

if ~isempty(res3)
    [VMAX,PHI] = meshgrid(linspace(min(phi3),max(phi3),200), linspace(min(vmax3),max(vmax3),200));
    subplot1(3)
    contourf(VMAX,PHI,griddata(phi3,vmax3,rmse3,VMAX,PHI),v)
    hold on
    for i=1:numel(ytic)-2,plot([xtic(1),xtic(end)],[ytic(i+1),ytic(i+1)],'-.','Color','black');end
    for i=1:numel(xtic)-2,plot([xtic(i+1),xtic(i+1)],[ytic(1),ytic(end)],'-.','Color','black');end
    plot([xtic(1),xtic(end),xtic(end),xtic(1),xtic(1)],[ytic(1),ytic(1),ytic(end),ytic(end),ytic(1)],...
    '-','LineWidth',2,'color','black')
    hold off
    set(gca,'fontsize', 12)
    caxis([min(v) max(v)])
    xticks(xtic)
    yticks(ytic)
    text(0.98,0.90,'SRTM3','FontSize',15,'FontWeight','bold',...
        'Units','normalized','HorizontalAlignment','right')
    xlabel('PHI [\circ]','FontWeight','bold')
    axis equal
end

cb = colorbar('south');
cb.Ticks = flipud(v');






% %%%%plots with vmax
% 
% figure
% dim = [5 10];%dimension of subplot
% subplot1(dim(1),dim(2),'Gap',[0.01 0.01],'XTickL','Margin','YTickL','Margin','FontS',10)
% for i=1:numel(phi_unique)
%     %extract needed data for subplot
%     idx_phi = find(phi == phi_unique(i));
%     vmax_plot = vmax(idx_phi);
%     rmse_plot = rmse(idx_phi);
%     %sort values according to phi
%     [vmax_plot,idx_sort] = sort(vmax_plot);
%     rmse_plot = rmse_plot(idx_sort);
%     %plot subplot
%     subplot1(i)
%     plot(vmax_plot,rmse_plot,'.-','MarkerSize',15)
%     [~,idx_min] = min(rmse_plot);
%     plot(vmax_plot(idx_min),rmse_plot(idx_min),'.','MarkerSize',15,'color','blue')
%     grid on
%     text(0.5,0.9,['PHI = ' num2str(phi_unique(i))],'FontSize',10,'FontWeight','bold',...
%         'Units','normalized','HorizontalAlignment','center')
%     if floor((i-1)/dim(2))+1==dim(1), xlabel('PHI','FontWeight','bold'); end
%     if mod(i-1,dim(2))==0, ylabel('RMSE','FontWeight','bold'); end
%     xlim([min(vmax_unique) max(vmax_unique)])
%     %ylim([0 500])
% end  

%%%%plots of vmax and phi together


end