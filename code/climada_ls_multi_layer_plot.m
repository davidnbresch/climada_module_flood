function climada_ls_multi_layer_plot(hazard,centroids,susc_check)

global climada_global
if ~climada_init_vars,  return;     end

if ~exist('hazard',     'var'), hazard      = [];   end
if ~exist('centroids',  'var'), centroids   = [];   end
if ~exist('susc_check', 'var'), susc_check  = 0;    end

if isempty(hazard)
    % prompt for RF hazard event set if not given
    hazard_set_file=[module_data_dir filesep 'hazards' filesep '*.mat'];
    [filename, pathname] = uigetfile(hazard_set_file,...
        'Select a land slide hazard event set:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        hazard_set_file =fullfile(pathname,filename);
    end
    load(hazard_set_file);
end

% max event
event_sum=sum(abs(hazard.intensity),2);
[~,sorted_i]=sort(event_sum,'descend');
event_i=sorted_i(1);

if any(centroids.lon~=hazard.lon)
    cprintf([1 0.5 0],'WARNING: centroids are incompatible with hazard, reconstructing centroids from hazard coords\n')
    centroids = [];
end

if isempty(centroids)
    centroids.lon           = hazard.lon;
    centroids.lat           = hazard.lat;
    centroids.centroid_ID   = hazard.centroid_ID;
    if isfield(hazard,'elevation_m')
        centroids.elevation_m = hazard.elevation_m;
    else
        centroids = climada_90m_DEM(centroids,'DL','NO_SAVE',4,0);
    end
    centroids               = centroids_TWI(centroids,0);
end

climada_figuresize(1,2)
view(60,30)

hold on

lon     = centroids.lon;
lat     = centroids.lat;
[x,y]   = meshgrid(unique(lon),unique(lat));

% mesh and normalise
tmp_TWI     = centroids.TWI;
tmp_TWI(~centroids.onLand) = NaN;
gridded_TWI = griddata(lon,lat,tmp_TWI,x,y,'linear');
gridded_TWI = gridded_TWI - min(min(gridded_TWI));
gridded_TWI = gridded_TWI./ max(max(gridded_TWI));
s(1) = surf(x,y,ones(size(gridded_TWI)).*4.5, 'edgecolor','none');
set(s(1),'CData',gridded_TWI,'FaceAlpha',0.9)
colormap(climada_colormap('FL',20))
freezeColors

tmp_slp     = centroids.slope_deg;
tmp_slp(~centroids.onLand) = NaN;
gridded_slp = griddata(lon,lat,tmp_slp,x,y,'linear');
gridded_slp = gridded_slp - min(min(gridded_slp));
gridded_slp = gridded_slp./ max(max(gridded_slp));
s(2) = surf(x,y,ones(size(gridded_slp)).*2.5, 'edgecolor','none');
set(s(2),'CData',gridded_slp,'FaceAlpha',0.9)
colormap(climada_colormap('schematic',20))
freezeColors

gridded_H2O = griddata(lon,lat,centroids.onLand,x,y,'nearest');
gridded_H2O(gridded_H2O==1) = NaN;
w = surf(x,y,gridded_H2O+0.05, 'edgecolor','none'); 
set(w,'FaceColor',getelements(seacolor(10),6),'FaceAlpha',0.9)
freezeColors

gridded_DEM = griddata(lon,lat,centroids.elevation_m,x,y,'linear');
gridded_DEM = gridded_DEM - min(min(gridded_DEM));
gridded_DEM = (gridded_DEM./ max(max(gridded_DEM))).*1; % factor for little gap between layers
s(3) = surf(x,y,gridded_DEM, 'edgecolor','none'); 
% colormap(getelements(landcolor(200),1:100))
colormap(gray)
[C, h] = contour3(x,y,gridded_DEM,5);
% set(s(1),'edgecolor','none');
l_h     = clabel(C,h);
for i=1:length(l_h)
    c_lbl = get(l_h(i),'String');
    c_lbl = str2num(c_lbl).*max(centroids.elevation_m) + min(centroids.elevation_m);
    c_lbl = sprintf('%4.0f',c_lbl);
    set(l_h(i),'String',c_lbl,'FontAngle','normal','Rotation',0,'FontWeight','bold');
end
for i =1:length(h), set(h(i),'EdgeColor','k','FaceColor','flat'); end
freezeColors

if susc_check
    if isfield(hazard,'IF_at_centroid')
        gridded_IF = griddata(lon,lat,hazard.IF_at_centroid,x,y,'nearest');
        gridded_DEM(gridded_IF==0) = NaN;
        s(4) = surf(x,y,gridded_DEM+0.05, 'edgecolor','none'); 
        set(s(4),'CData',gridded_IF,'FaceAlpha',1)
        layer_4_name = 'Relative susceptibility';
        caxis([min(hazard.IF_at_centroid) max(hazard.IF_at_centroid)])
    else
        hazard.factor_of_safety(:,~centroids.onLand) = NaN;
        FoS = min(hazard.factor_of_safety,[],1);
    %     gridded_FoS = griddata(lon,lat,hazard.factor_of_safety(event_i,:),x,y,'linear');
        gridded_FoS = griddata(lon,lat,FoS,x,y,'linear');
    %     gridded_FoS(gridded_FoS > 1 | gridded_FoS < 0) = NaN;
        s(4) = surf(x,y,ones(size(gridded_FoS)).*3, 'edgecolor','none');
        set(s(4),'CData',gridded_FoS,'FaceAlpha',0.6)
        layer_4_name = 'Factor of safety';
    end
    colormap(getelements(flipud(hot(30)),10:30))
    cbar = colorbar;
    set(get(cbar,'ylabel'),'String', {layer_4_name} ,'fontsize',10);
    set(cbar,'YTickLabel',[]);
    freezeColors
end


% l = legend(s(1:3),'Topographic wetness index','Steepness','Digital elevation model');
% set(l,'box','off','Location','NorthEast')
% shading interp

ylabel('Latitude')
xlabel('Longitude')
set(gca,'ZTick',[], 'ZColor', 'w');
axis tight
axis([min(lon) max(lon) min(lat) max(lat)])
axis off
title({'Multi-layered Approach to Land Slide Simulation',''},'fontsize',18)


