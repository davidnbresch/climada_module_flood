function script_FPcalibration
%%%%toDo
%in this script: include calibration for all 3 DEMs
%test climada_calibration_calcRMSE (with rmv-vector) and include in this script
%finish restructering of climada_calibration_ls (include
%savepath,documentation...)
%rewrite calibration_snap such that it does work without startarea_endarea

%%
%define original data and save directories
dat_dir = 'C:\Users\Simon Rölli\Desktop\data\calibration';

%dem which shall be calibration high resolution (org) and low resolution (snap)
file_dem_org = 'C:\Users\Simon Rölli\Desktop\data\centroids_large\_LS_Sarnen_alti3d_2m.mat';
org_centroids = load(file_dem_org,'centroids');
org_centroids = org_centroids.centroids;

orgS = shaperead('C:\Users\Simon Rölli\Desktop\data\inventory\be_polygon_forMatlab.shp');
subS = shaperead('C:\Users\Simon Rölli\Desktop\data\calibration\subData\subS_2x3m.shp');

%%
%load data DEM and snap for all DEMs

%srtm3
file_dem_snap = 'C:\Users\Simon Rölli\Desktop\data\centroids_large\_LS_Sarnen_srtm3'; %srtm3
cal_centroids_srtm3 = load(file_dem_snap,'centroids');
cal_centroids_srtm3 = cal_centroids_srtm3.centroids;
snapS_srtm3 = shaperead('C:\Users\Simon Rölli\Desktop\data\calibration\snapData\snapS_63x92m.shp');

%srtm1
file_dem_snap = 'C:\Users\Simon Rölli\Desktop\data\centroids_large\_LS_Sarnen_srtm1';%srtm1
cal_centroids_srtm1 = load(file_dem_snap,'centroids');
cal_centroids_srtm1 = cal_centroids_srtm1.centroids;
snapS_srtm1 = shaperead('C:\Users\Simon Rölli\Desktop\data\calibration\snapData\snapS_21x30m.shp');

%alti3d
file_dem_snap = 'C:\Users\Simon Rölli\Desktop\data\centroids_large\_LS_Sarnen_interpol_alti3d';%alti3d ca. 7*10m
cal_centroids_alti3d = load(file_dem_snap,'centroids');
cal_centroids_alti3d = cal_centroids_alti3d.centroids;
snapS_alti3d = shaperead('C:\Users\Simon Rölli\Desktop\data\calibration\snapData\snapS_7x10m.shp');

%%
%create parameter sets
dH = 0;
phi = [10:50];
vmax = [1:12];
exp = 25;
iT = 0.0003;
perWT = [1 0.8 0.4 0 0 0 0.4 0.8];

%create all combinations
sets = {dH,vmax,phi,exp,iT};
[dH vmax phi exp iT] = ndgrid(sets{:});
dH = dH(:); phi = phi(:); vmax = vmax(:); exp = exp(:); iT = iT(:);
perWT = repmat(perWT,numel(phi),1);

flowPara.dH = dH;
flowPara.exp = exp;
flowPara.phi = phi;
flowPara.vmax = vmax;
flowPara.iT = iT;
flowPara.perWT = perWT;
%%
%run model for calibration

[res_files_srtm3,unit_srtm3] = climada_calibration_ls(org_centroids,cal_centroids_srtm3,orgS,flowPara,dat_dir,subS,snapS_srtm3);
[res_files_srtm1,unit_srtm1]= climada_calibration_ls(org_centroids,cal_centroids_srtm1,orgS,flowPara,dat_dir,subS,snapS_srtm1);
[res_files_alti3d,unit_alti3d] = climada_calibration_ls(org_centroids,cal_centroids_alti3d,orgS,flowPara,dat_dir,subS,snapS_alti3d);

%%
%select slide to remove

%set some thresholds
max_area = 30000; %[m^2] according to Millage etal area (10^1–10^4 m2).
max_lgt = 1000;

%extract data from subS
lgt = [subS.length];
area = [subS.area];

rmv_srtm3 = logical([snapS_srtm3.removed]);
rmv_srtm1 = logical([snapS_srtm1.removed]);
rmv_alti3d = logical([snapS_alti3d.removed]);

rmv_srtm3(lgt>max_lgt) = 1;
rmv_srtm1(lgt>max_lgt) = 1;
rmv_alti3d(lgt>max_lgt) = 1;

rmv_srtm3(area>max_area) = 1;
rmv_srtm1(area>max_area) = 1;
rmv_alti3d(area>max_area) = 1;



%%
if 0
    %calculate RMSE full 
    rmse_srtm3 = climada_calibration_calcRMSE(res_files_srtm3,subS);
    rmse_srtm1 = climada_calibration_calcRMSE(res_files_srtm1,subS);
    rmse_alti3d = climada_calibration_calcRMSE(res_files_alti3d,subS);

    disp('hier')
    %calculate RMSE excluding removed slides
    rmse_srtm3_rmv = climada_calibration_calcRMSE(res_files_srtm3,subS,[],rmv_srtm3);
    rmse_srtm1_rmv = climada_calibration_calcRMSE(res_files_srtm1,subS,[],rmv_srtm1);
    rmse_alti3d_rmv = climada_calibration_calcRMSE(res_files_alti3d,subS,[],rmv_alti3d);
    
    %calculate RMSE excluding slides with lower local max slope than phi
    rmse_srtm3_maxslope = climada_calibration_calcRMSE(res_files_srtm3,subS,snapS_srtm3);
    rmse_srtm1_maxslope = climada_calibration_calcRMSE(res_files_srtm1,subS,snapS_srtm1);
    rmse_alti3d_maxslope = climada_calibration_calcRMSE(res_files_alti3d,subS,snapS_alti3d);
    
    %calculate RMSE excluding slides with lower local max slope than phi
    %and exclude removed
    rmse_srtm3_both = climada_calibration_calcRMSE(res_files_srtm3,subS,snapS_srtm3,rmv_srtm3);
    rmse_srtm1_both = climada_calibration_calcRMSE(res_files_srtm1,subS,snapS_srtm1,rmv_srtm1);
    rmse_alti3d_both = climada_calibration_calcRMSE(res_files_alti3d,subS,snapS_alti3d,rmv_alti3d);
    
    rmv_srtm3(area<unit_srtm3) = 1;
    rmv_srtm1(area<unit_srtm1) = 1;
    rmv_alti3d(area<unit_alti3d) = 1;
    
    %calculate RMSE excluding removed slides (additionally removed small
    %slides --> smaller then unitarea)
    rmse_srtm3_rmvsmall = climada_calibration_calcRMSE(res_files_srtm3,subS,[],rmv_srtm3);
    rmse_srtm1_rmvsmall = climada_calibration_calcRMSE(res_files_srtm1,subS,[],rmv_srtm1);
    rmse_alti3d_rmvsmall = climada_calibration_calcRMSE(res_files_alti3d,subS,[],rmv_alti3d);
    
    %calculate RMSE excluding everything: removed slides (additionally removed small
    %slides --> smaller then unitarea) and slides with local slope smaller
    %than phi
    rmse_srtm3_bothsmall = climada_calibration_calcRMSE(res_files_srtm3,subS,snapS_srtm3,rmv_srtm3);
    rmse_srtm1_bothsmall = climada_calibration_calcRMSE(res_files_srtm1,subS,snapS_srtm1,rmv_srtm1);
    rmse_alti3d_bothsmall = climada_calibration_calcRMSE(res_files_alti3d,subS,snapS_alti3d,rmv_alti3d);
    clear res_files_srtm3 res_files_srtm1 res_files_alti3d
    
    save([dat_dir filesep 'rmse_results.mat'],'rmse_srtm3','rmse_srtm1','rmse_alti3d',...
    'rmse_srtm3_rmv','rmse_srtm1_rmv','rmse_alti3d_rmv','rmse_srtm3_maxslope',...
    'rmse_srtm1_maxslope','rmse_alti3d_maxslope','rmse_srtm3_both','rmse_srtm1_both',...
    'rmse_alti3d_both','rmse_srtm3_rmvsmall','rmse_srtm1_rmvsmall','rmse_alti3d_rmvsmall',...
    'rmse_srtm3_bothsmall','rmse_srtm1_bothsmall','rmse_alti3d_bothsmall')
end
%%
%plot RMSE
load([dat_dir filesep 'rmse_results.mat'],'rmse_srtm3','rmse_srtm1','rmse_alti3d',...
    'rmse_srtm3_rmv','rmse_srtm1_rmv','rmse_alti3d_rmv','rmse_srtm3_maxslope',...
    'rmse_srtm1_maxslope','rmse_alti3d_maxslope','rmse_srtm3_both','rmse_srtm1_both',...
    'rmse_alti3d_both','rmse_srtm3_rmvsmall','rmse_srtm1_rmvsmall','rmse_alti3d_rmvsmall',...
    'rmse_srtm3_bothsmall','rmse_srtm1_bothsmall','rmse_alti3d_bothsmall')


% %plot all vs rmv
% climada_calibration_plotRMSE(rmse_alti3d,rmse_srtm1,rmse_srtm3,rmse_alti3d_rmv,rmse_srtm1_rmv,rmse_srtm3_rmv);
% %plot all vs maxslope
% climada_calibration_plotRMSE(rmse_alti3d,rmse_srtm1,rmse_srtm3,rmse_alti3d_maxslope,rmse_srtm1_maxslope,rmse_srtm3_maxslope);
% %plot all vs both
% climada_calibration_plotRMSE(rmse_alti3d,rmse_srtm1,rmse_srtm3,rmse_alti3d_both,rmse_srtm1_both,rmse_srtm3_both);
% %plot all vs removed small 
% climada_calibration_plotRMSE(rmse_alti3d,rmse_srtm1,rmse_srtm3,rmse_alti3d_rmvsmall,rmse_srtm1_rmvsmall,rmse_srtm3_rmvsmall);
% %plot all vs both_removed small 
% climada_calibration_plotRMSE(rmse_alti3d,rmse_srtm1,rmse_srtm3,rmse_alti3d_bothsmall,rmse_srtm1_bothsmall,rmse_srtm3_bothsmall);

%plots for thesis
climada_calibration_plotRMSE(rmse_alti3d_rmv,rmse_srtm1_rmv,rmse_srtm3_rmv)

%plots for appendix
climada_calibration_plotRMSE(rmse_alti3d_maxslope,rmse_srtm1_maxslope,rmse_srtm3_maxslope)

%plot of number of slides considered --> only for maxslope possible (other
%with rmv slides --> number of observations is not taken
% rmse_srtm3_maxslope_num = rmse_srtm3_maxslope;
% c = num2cell([rmse_srtm3_maxslope_num.num_obs]);
% [rmse_srtm3_maxslope_num.rmse_lgt] = c{:};
% rmse_srtm1_maxslope_num = rmse_srtm1_maxslope;
% c = num2cell([rmse_srtm1_maxslope_num.num_obs]);
% [rmse_srtm1_maxslope_num.rmse_lgt] = c{:};
% rmse_alti3d_maxslope_num = rmse_alti3d_maxslope;
% c = num2cell([rmse_alti3d_maxslope_num.num_obs]);
% [rmse_alti3d_maxslope_num.rmse_lgt] = c{:};
% climada_calibration_plotRMSE(rmse_alti3d_maxslope,rmse_srtm1_maxslope,rmse_srtm3_maxslope,rmse_alti3d_maxslope_num,rmse_srtm1_maxslope_num,rmse_srtm3_maxslope_num)

%remove small
climada_calibration_plotRMSE(rmse_alti3d_rmvsmall,rmse_srtm1_rmvsmall,rmse_srtm3_rmvsmall)
%climada_calibration_plotRMSE(rmse_alti3d_bothsmall,rmse_srtm1_bothsmall,rmse_srtm3_bothsmall)
%area --> code needed to be changed
%[f1,f2]=climada_calibration_plotRMSE(rmse_alti3d_rmv,rmse_srtm1_rmv,rmse_srtm3_rmv)



%plot which shows distribution of source cell slope
% subS = shaperead('C:\Users\Simon Rölli\Desktop\data\calibration\subData\subS_2x3m.shp');
% snapS_srtm3 = shaperead('C:\Users\Simon Rölli\Desktop\data\calibration\snapData\snapS_63x92m.shp');
% snapS_srtm1 = shaperead('C:\Users\Simon Rölli\Desktop\data\calibration\snapData\snapS_21x30m.shp');
% snapS_alti3d = shaperead('C:\Users\Simon Rölli\Desktop\data\calibration\snapData\snapS_7x10m.shp');
% 
% figure
% b=boxplot([[subS.slope];[snapS_alti3d.slope];[snapS_srtm1.slope];[snapS_srtm3.slope]]',...
%     {'ALTI3D (2m)','ALTI3D (10m)','SRTM1 (30m)','SRTM3 (90m)'},'Widths',[0.2 0.2 0.2 0.2],...
%     'FactorGap',0.001)
% ylabel('Slope of Source Cell [\circ]')
% yl = ylim;
% ylim([0 yl(2)])
% xl = xlim; yl = ylim;
% hold on
% plot([xl(1) xl(2) xl(2) xl(1) xl(1)],[yl(1) yl(1) yl(2) yl(2) yl(1)],'-','LineWidth',2,'Color','black')
% 
% hold on

%plot which shows xy plot (modelled lgt vs observed lgt) for phi =27, vmax
%=1

subS = shaperead('C:\Users\Simon Rölli\Desktop\data\calibration\subData\subS_2x3m.shp');
snapS_srtm3 = shaperead('C:\Users\Simon Rölli\Desktop\data\calibration\snapData\snapS_63x92m.shp');
snapS_srtm1 = shaperead('C:\Users\Simon Rölli\Desktop\data\calibration\snapData\snapS_21x30m.shp');
snapS_alti3d = shaperead('C:\Users\Simon Rölli\Desktop\data\calibration\snapData\snapS_7x10m.shp');
PHI27_alti3d = shaperead('C:\Users\Simon Rölli\Desktop\data\calibration\caliData\7x10m\7x10m_phi27_vmax4_exp25_dH0_iT00003_perWT1_08_04_0_0_0_04_08.shp');
PHI27_srtm1 = shaperead('C:\Users\Simon Rölli\Desktop\data\calibration\caliData\21x30m\21x30m_phi27_vmax4_exp25_dH0_iT00003_perWT1_08_04_0_0_0_04_08.shp');
PHI27_srtm3 = shaperead('C:\Users\Simon Rölli\Desktop\data\calibration\caliData\63x92m\63x92m_phi27_vmax4_exp25_dH0_iT00003_perWT1_08_04_0_0_0_04_08.shp');

figure('units','normalized','outerposition',[0 0.2 0.7 0.7])
subplot1(1,2,'Gap',[0.03 0.01],'FontS',12)

%read out data
obs_lgt = [subS.length];
mod_lgt = [PHI27_alti3d.length];
max_srcslope = [snapS_alti3d.max_srcslop];
idx = find(max_srcslope<27);

subplot1(1)
plot(obs_lgt,mod_lgt,'.','MarkerSize',10)
grid on
ylabel('modelled length [m] with PHI=27\circ and VMAX=4m/s','FontWeight','bold')
xlabel('observed length [m]','FontWeight','bold')
yl = xlim;
ylim(yl);
xt = yticks;
xticks(xt);
hold on
plot(obs_lgt(idx),mod_lgt(idx),'.','color','red','MarkerSize',10)
plot([yl(1) yl(2)],[yl(1) yl(2)],'--','color','black')
hold off
subplot1(2)
plot(max_srcslope,mod_lgt,'.','MarkerSize',10)
grid on
xlabel('maximum slope gradient of source cell [\circ]','FontWeight','bold')
ylim(yl);
hold on
plot(max_srcslope(idx),mod_lgt(idx),'.','color','red','MarkerSize',10)
hold off






end