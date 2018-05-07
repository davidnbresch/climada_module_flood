function [S,spreaded,IDfield] = climada_ls_FPcalibration(orgS,orgLon,orgLat,orgStart,orgEnd,...
    lon,lat,elevation,field)

% Script to calibrate flow path parameters
% MODULE:
%   flood
% NAME:
%   climada_ls_FPcalibration
% PURPOSE:
%   
% CALLING SEQUENCE:
%   
% EXAMPLE:
%   
% INPUTS: 
%     calculate: 1 if area lenght start end need to be calcualted
%     calculate2: 1 if 
% OPTIONAL INPUT PARAMETERS:
%    
% OUTPUTS:
%    
% MODIFICATION HISTORY:
% Thomas Rölli, thomasroelli@gmail.com, 20180502, init

global climada_global
if ~climada_init_vars, return; end

if ~exist('calculate') calculate=0; end
if ~exist('calculate2') calculate2=0; end

%init S
S = [];

n_lon = size(lon,2);
n_lat = size(lat,1);
start = zeros(n_lat,n_lon);
end_ = zeros(n_lat,n_lon);

%caculation of unitarea --> not considering slope
deg_km = 111.32;
dlat = abs(min(diff(lat(:,1)))); 
dlon = abs(min(diff(lon(1,:))));
dy = dlat*(deg_km * 1000);
dx = dlon*cosd(mean(lat(:,1)))*(deg_km * 1000); 
unitarea = dx*dy;

%prepare grid to use in dsearchn()
grid = [lon(:) lat(:)];

%when get to coarser resolution --> 
%indices of start areas
fprintf('Find start cells of new/coarser grid\n');
idx_st = find(orgStart~=0); %indices of original start
points = [orgLon(idx_st) orgLat(idx_st)];
k_st = dsearchn(grid,points); 

fprintf('Find end cells of new/coarser grid\n');
idx_en = find(orgEnd~=0); %indices of original end
points = [orgLon(idx_en) orgLat(idx_en)]; %original points
k_en = dsearchn(grid,points); %get indices of nearest point to original points

maxfield='AREA_GIS';
start = zeros(size(elevation));
end_ = zeros(size(elevation));
length = zeros(size(elevation)); %to check if other (longer) slide exists at same cell already
for i=1:numel(orgS)
    S(i).(field) = orgS(i).(field); %transform field ID
    
    ik_st = find(orgStart(idx_st) == orgS(i).(field)); %find corresponding index of field ID
    ik_en = find(orgEnd(idx_en) == orgS(i).(field));
    %write start and end coordinates in Structure
    S(i).startlon = lon(k_st(ik_st));
    S(i).startlat = lat(k_st(ik_st));
    S(i).endlon = lon(k_en(ik_en));
    S(i).endlat = lat(k_en(ik_en));
    
    if i==65 
        disp('')
    end
    %calculate distance of slide in new raster (with slope)
    [Ist,Jst] = ind2sub(size(lon),k_st(ik_st));
    [Ien,Jen] = ind2sub(size(lon),k_en(ik_en));
    dI = abs(Ist-Ien)*dy;
    dJ = abs(Jst-Jen)*dx;
    dz = elevation(k_st(ik_st))-elevation(k_en(ik_en));
    try 
        lgt = double(sqrt(sqrt(dI^2+dJ^2)^2+dz^2));
    catch
        lgt = double(0);
    end
    %write snapped length in structure
    S(i).snap_length = lgt;
    
    %create start matrix --> if one cell is affected several times-->
    %take the one which is longer (of original raster length)
    %for end matrix not so important --> just overwrite
    old_lgt = length(k_st(ik_st));
    if old_lgt < orgS(i).R_LGT_SL
        length(k_st(ik_st)) = orgS(i).R_LGT_SL;
        start(k_st(ik_st)) = orgS(i).ID;
    end
    end_(k_en(ik_en)) = orgS(i).ID;
    
end 


%area of each raster cell --> considering slope
cell_area = climada_centroids_area(lon,lat,elevation,0);

%parameters for flowpath
dH = 0;
exponent = 25;
v_max = 4;
phi = 22;
delta_i = 0.0003;
perWt = [1 0.8 0.4 0 0 0 0.4 0.8];

%needed data for flowpath
source = zeros(size(elevation));
elevation = fillsinks(elevation);
mult_flow = climada_ls_multipleflow(lon,lat,elevation,exponent,dH,1);
[~,horDist,verDist] = climada_centroids_gradients(lon,lat,elevation);

%init of matrix to save flow paths ID and log(maximum intensity)
slides_ID = zeros(size(elevation));
slides_int = zeros(size(elevation));



source = zeros(size(elevation));
source(find(start~=0)) = 1;
sel_source = zeros(size(elevation));
mask = zeros(size(elevation));
spreaded = zeros(size(elevation)); %saves intensity (max when several slides)
temp_spreaded = zeros(size(elevation));
IDfield = zeros(size(elevation))+100000; %saves ID of slides (sum when several slides)

lgt = zeros(size(S));
area = zeros(size(S));



k=0;

%%
%method which chooses several slide sources at a time to save computing
%time --> is done by taking slides which are far enough apart from each
%other (distance defined by buf_m)
numsource = find(source==1);
msgstr   = sprintf('Assessing flow of %i flat cells ... ',numel(numsource));
mod_step = 10; % first time estimate after 10 assets, then every 100
if climada_global.waitbar
    fprintf('%s (updating waitbar with estimation of time remaining every 100th centroid)\n',msgstr);
    h        = waitbar(0,msgstr);
    set(h,'Name','Assigning TWI');
else
    fprintf('%s (waitbar suppressed)\n',msgstr);
    format_str='%s';
end

buf_m = 1000;
imask = ceil(buf_m/dy);
jmask = ceil(buf_m/dx);
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

    %spread selected source cells
    temp_spreaded = climada_ls_propagation(sel_source,mult_flow,horDist,verDist,v_max,phi,delta_i,perWt);
    
    %%
    %reconnect spreaded slides and derive area and lgt
    temp_source_ID = zeros(size(elevation));
    temp_source_ID = sel_source.*start;
    [temp_IDfield,temp_ar,temp_lgt,tempID] = climada_ls_slideScores(temp_spreaded,temp_source_ID,elevation,cell_area,dx,dy);
    
    temp_IDfield(temp_IDfield==0) = 100000;
    
    %find corresponding index in S of slides in sel_source; considering
    %'field' and save in array
    c = bsxfun(@eq,tempID,[S.(field)]');
    [sidx,~] = find(c); %gives field number (ID) of slide which shall be use in S.(field) to assign values
    %write lenght and area in vector
    lgt(sidx) = temp_lgt;
    area(sidx) = temp_ar;
    %merge intensity (max when overlapping) and ID (sum...)
    spreaded = max(spreaded,temp_spreaded);
    IDfield = min(IDfield,temp_IDfield);
   
    %set to zero for next round of selection
    mask(:) = 0;
    sel_source(:) = 0;
    
    %progress
     if mod(i,mod_step)==0
            mod_step = 1000;
            msgstr = sprintf('%i/%i spreaded',numel(numsource-numel(find(source==1)),numsource));
            if climada_global.waitbar
                waitbar(i/numel(S),h,msgstr); % update waitbar
            else
                fprintf(format_str,msgstr); % write progress to stdout
                format_str=[repmat('\b',1,length(msgstr)) '%s']; % back to begin of line
            end
        end
end %while --> ends when no sources left

c = num2cell(area);
[S.trig_area] = c{:};

c = num2cell(lgt);
[S.trig_lgt] = c{:};

IDfield(IDfield==100000) = 0;


end