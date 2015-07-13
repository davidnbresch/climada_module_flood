function centroids=centroids_LAI(centroids, check_plots)
% Assign Leaf Area Indices (LAIs) to given centroids
% MODULE:
%   flood
% NAME:
%	centroids_LAI_assign
% PURPOSE:
%   Determine Leaf Area Indices (LAIs) for given centroids by
%   reading a grayscale picture of global annual mean LAIs.
%   Leaf area index (LAI) is a dimensionless quantity that characterizes
%   plant canopies. It is defined as the one-sided green leaf area per unit
%   ground surface area (LAI = leaf area / ground area, m^2/m^2).
%   For more information on the LAI see
%       http://en.wikipedia.org/wiki/Leaf_area_index
%   The Moderate Resolution Imaging Spectroradiometer (MODIS) aboard NASA's
%   Terra and Aqua satellites collects global LAI data on a daily basis.
%   Values range from 0 to 7 square meters of leaf area per square meter of
%   land surface.
%   Global monthly mean LAI images in grayscale can be downloaded here:
%       http://neo.sci.gsfc.nasa.gov/view.php?datasetId=MOD15A2_M_LAI
% NOTE: Since only monthly mean data are provided under the link above, a
% grayscale picture of the annual mean LAI needs to be produced manually by
% adding monthly pictures and taking the average.
%
% CALLING SEQUENCE:
%   centroids = centroids_LAI(centroids,check_plots)
% EXAMPLE:
%   centroids = centroids_LAI(centroids)
% INPUTS:
%  centroids: Climada centroids struct; the following fields are required:
%         .lat:           Latitude
%         .lon:           Longitude
%         .centroid_ID:   centroid ID
%         .admin0_name    Country name
% OPTIONAL INPUT PARAMETERS:
%   LAI_img_filename: File name of the global mean LAI image. If empty, a
%   default image is used, and if it does not exists, the function prompts
%   for the image filename.
%   check_plots: whether a plot should be drawn (=1), for each month (=2) or not at all (=0; default) 
%                   [about 3x slower when using check_plots = 2]
% OUTPUTS:
%   centroids: The input centroids structure with an additional field
%   'leaf_area_index', which contains the Leaf Area Index for each centroid
% NOTE:
% MODIFICATION HISTORY:
% Melanie Bieli, melanie.bieli@bluewin.ch, 20150318, initial
% Gilles Stassen, gillesstassen@hotmail.com, 20150407, monthly data; auto
%                   file finder/downloader; interp2 to determine values at centroids; plotter
%-


% import/setup global variables
global climada_global
if ~climada_init_vars,return;end;

% check arguments
if ~exist('centroids','var') || isempty(centroids)
    climada_centroids_load
end
if ~exist('LAI_img_filename','var'),LAI_img_filename=''; end
if ~exist('check_plots', 'var') || isempty(check_plots)
    check_plots = 0;    end


% locate the module's data folder
module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];
% http://neo.sci.gsfc.nasa.gov/servlet/RenderData?si=1610438&cs=gs&format=TIFF&width=3600&height=1800
% ftp://neoftp.sci.gsfc.nasa.gov/gs/MOD15A2_M_LAI/MOD15A2_M_LAI_2014-01.TIFF
% PARAMETERS
%
% the file with the LAI data

min_lat=-90; % degree, defined on the webpage above
max_lat= 90; % defined on the webpage above
%
% prepare bounding box 
bbox = [min(centroids.lon)-1, min(centroids.lat)-1,max(centroids.lon)+1, max(centroids.lat)+1];
% for progress mgmt
t0 = clock;
format_str = '%s';

if check_plots
    fig = figure('name','Leaf area index','Color','w');
end

for month_i = 1 : 12
    LAI_img_fN = sprintf('MOD15A2_M_LAI_%04d-%02d.PNG',str2num(datestr(now,'yyyy'))-1,month_i);
    LAI_img_URL = ['ftp://neoftp.sci.gsfc.nasa.gov/gs/MOD15A2_M_LAI/' LAI_img_fN];
    LAI_img_fullpath =[module_data_dir filesep 'system' filesep 'MOD15A2_M_LAI' filesep LAI_img_fN];
    
    [fP,fN,~] = fileparts(LAI_img_fullpath);
    
    if ~exist(fP,'dir')
        mkdir(fP);
    end
    
    LAI_img_fullpath_mat = [fP filesep fN '.mat'];
    
    if climada_check_matfile(LAI_img_fullpath_mat) % try matfile
        load(LAI_img_fullpath_mat);
    elseif exist(LAI_img_fullpath,'file') % try original png file
        LAI_img_filename = LAI_img_fullpath;
        full_img=imread(LAI_img_fullpath);
        save(LAI_img_fullpath_mat,'full_img');
    else % look for png file using subdir
        LAI_img_dir = subdir([climada_global.root_dir filesep LAI_img_fN]);
        if ~isempty(LAI_img_dir)
            s = movefile(LAI_img_dir.name,fP);
            full_img=imread(LAI_img_fullpath);
            save(LAI_img_fullpath_mat,'full_img')
        else % not found using subdir, read from URL
            full_img=imread(LAI_img_URL);
            save(LAI_img_fullpath_mat,'full_img');
        end
    end
    
    % White means no data available, hence we set all white pixels to 0
    % (i.e., to black)
    full_img(full_img==255) = NaN;
    full_img = flipud(full_img);
    
    [X, Y] = meshgrid([-180:0.1:179.9],[-90:0.1:89.9]);
    
    if (check_plots == 1 && month_i == 7) || check_plots ==2
        figure(fig)
        colormap(LAI_colormap);
        imagesc([-180:0.1:179.9],[-90:0.1:89.9],full_img.*(7/255));
        hold on
        set(gca,'ydir','normal');
        shading flat
        climada_plot_world_borders
        title_str = sprintf('Leaf Area Index (m^2/m^2) %s',...
            datestr(datenum(sprintf('%02d',month_i),'mm'),'mmm'));
        if isfield(centroids,'admin0_name')
            title_str = [title_str ' - ' centroids.admin0_name];
        end
        title(title_str)
        colorbar
        % axis([min(centroids.lon) max(centroids.lon) min(centroids.lat) max(centroids.lat)])
        hold off
        drawnow
    end
    
    centroids.LAI(month_i,:)=interp2(X, Y, double(full_img),centroids.lon,centroids.lat,'linear');
    
    t_elapsed_month   = etime(clock,t0)/month_i;
    months_remaining  = 12-month_i;
    t_projected_sec   = t_elapsed_month*months_remaining;
    if t_projected_sec<60
        msgstr = sprintf('est. %3.0f sec left (%i/%i months)',t_projected_sec,month_i,12);
    else
        msgstr = sprintf('est. %3.1f min left (%i/%i months)',t_projected_sec/60,month_i,12);
    end
    fprintf(format_str,msgstr); % write progress to stdout
    format_str=[repmat('\b',1,length(msgstr)) '%s']; % back to begin of line
end
fprintf(format_str, sprintf('processing leaf area index took %2.1f seconds \n',etime(clock,t0)));

% convert range of greyscale values (0 to 255) into actual Leaf Area
% Indices (LAIs) ranging from 0 to 7
centroids.LAI = centroids.LAI.*(7/255);

end

function map = LAI_colormap
map = [ 0.966413 0.987374 0.991665
        0.961984 0.985652 0.990681
        0.957555 0.983929 0.989696
        0.953126 0.982207 0.988712
        0.948697 0.980484 0.987728
        0.944268 0.978762 0.986744
        0.939839 0.977040 0.985759
        0.935409 0.975317 0.984775
        0.930980 0.973595 0.983791
        0.928766 0.972734 0.983299
        0.922122 0.970150 0.981822
        0.919908 0.969289 0.981330
        0.913264 0.966705 0.979854
        0.911050 0.965844 0.979362
        0.904406 0.963260 0.977885
        0.902191 0.962399 0.977393
        0.894579 0.959539 0.973841
        0.888428 0.957324 0.969166
        0.885352 0.956217 0.966828
        0.876125 0.952895 0.959815
        0.869973 0.950681 0.955140
        0.863822 0.948466 0.950465
        0.860746 0.947359 0.948128
        0.851519 0.944037 0.941115
        0.845367 0.941822 0.936440
        0.839216 0.939608 0.931765
        0.836140 0.938501 0.929427
        0.826913 0.935179 0.922414
        0.820761 0.932964 0.917739
        0.814610 0.930750 0.913064
        0.811534 0.929642 0.910727
        0.802307 0.926321 0.903714
        0.792157 0.922414 0.897501
        0.779608 0.917493 0.890365
        0.767059 0.912572 0.883230
        0.754510 0.907651 0.876094
        0.741961 0.902730 0.868958
        0.735686 0.900269 0.865390
        0.716863 0.892887 0.854687
        0.704314 0.887966 0.847551
        0.691765 0.883045 0.840415
        0.679216 0.878124 0.833280
        0.666667 0.873203 0.826144
        0.654118 0.868281 0.819008
        0.641569 0.863360 0.811872
        0.635294 0.860900 0.808305
        0.616471 0.853518 0.797601
        0.603922 0.848597 0.790465
        0.591373 0.843337 0.781976
        0.578824 0.837924 0.772872
        0.566275 0.832511 0.763768
        0.553726 0.827097 0.754664
        0.541176 0.821684 0.745559
        0.534902 0.818977 0.741007
        0.516078 0.810857 0.727351
        0.503529 0.805444 0.718247
        0.490980 0.800031 0.709143
        0.478431 0.794617 0.700038
        0.465882 0.789204 0.690934
        0.453333 0.783791 0.681830
        0.440784 0.778378 0.672726
        0.434510 0.775671 0.668174
        0.415686 0.767551 0.654518
        0.403137 0.762138 0.645413
        0.393172 0.757093 0.634648
        0.384068 0.752172 0.623330
        0.374963 0.747251 0.612011
        0.365859 0.742330 0.600692
        0.356755 0.737409 0.589373
        0.347651 0.732488 0.578055
        0.338547 0.727566 0.566736
        0.329443 0.722645 0.555417
        0.320338 0.717724 0.544098
        0.311234 0.712803 0.532780
        0.302130 0.707882 0.521461
        0.297578 0.705421 0.515802
        0.283922 0.698039 0.498824
        0.274817 0.693118 0.487505
        0.265713 0.688197 0.476186
        0.256609 0.683276 0.464867
        0.248904 0.675356 0.452949
        0.241523 0.666744 0.440892
        0.234141 0.658132 0.428835
        0.226759 0.649519 0.416778
        0.219377 0.640907 0.404721
        0.211995 0.632295 0.392664
        0.204614 0.623683 0.380607
        0.197232 0.615071 0.368551
        0.189850 0.606459 0.356494
        0.182468 0.597847 0.344437
        0.175087 0.589235 0.332380
        0.171396 0.584929 0.326351
        0.160323 0.572011 0.308266
        0.152941 0.563399 0.296209
        0.145559 0.554787 0.284152
        0.138178 0.546175 0.272095
        0.129719 0.538639 0.265206
        0.121107 0.531257 0.259054
        0.112495 0.523875 0.252903
        0.103883 0.516494 0.246751
        0.095271 0.509112 0.240600
        0.086659 0.501730 0.234448
        0.078047 0.494348 0.228297
        0.069435 0.486967 0.222145
        0.060823 0.479585 0.215994
        0.052211 0.472203 0.209842
        0.043599 0.464821 0.203691
        0.039293 0.461130 0.200615
        0.026374 0.450058 0.191388
        0.017762 0.442676 0.185236
        0.009150 0.435294 0.179085
        0.000538 0.427912 0.172933
        0.000000 0.417993 0.168627
        0.000000 0.407905 0.164444
        0.000000 0.397816 0.160261
        0.000000 0.387728 0.156078
        0.000000 0.377639 0.151895
        0.000000 0.367551 0.147712
        0.000000 0.357463 0.143529
        0.000000 0.347374 0.139346
        0.000000 0.337286 0.135163
        0.000000 0.327197 0.130980
        0.000000 0.317109 0.126797
        0.000000 0.312065 0.124706
        0.000000 0.296932 0.118431
        0.000000 0.286844 0.114248
        0.000000 0.276755 0.110065
        0.000000 0.266667 0.105882 ];
end