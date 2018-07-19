function climada_ls_probAssPlot(modS,obsS,fieldname,qqOrHist,limQQ,binWidth,limHIST,CDF_plot,limCDF)

% Plots distributions (qq-plot and histogram) of probabilistically modelled and
% observed (landslide inventory) variable 
% Important: structure of input data is too complicated at moment --> work
% in progress. --> at moment a lot of effort is needed to get data to the
% right structure
% qq --> several qq plot within subplot (for each dem) possible if several
% datasets are provided (e.g. modS.SRTM3.data1, modS.SRTM3.data2 and same
% for obsS.SRTM3.data1....)
% histogram --> just first data couple is plotted (e.g. modS.SRTM3.data1
% vs obsS.SRTM3.data1) for each DEM in a subplot
% cdf --> different --> just modS is plotted --> extract data from
% modS.dem.data1..data2  and plots as single cdf line; for obsS provide the
% same dataset as modS --> qq plot will be on 1:1 line
% MODULE:
%  
% NAME:
%  climada_probSet_ls
% PURPOSE:
%  Different plots which describe distribution: qq-plot --> plots each modS
%  and obsS pair for all DEM
% CALLING SEQUENCE:
%   
% EXAMPLE:
%   
% INPUTS:
%   modS:   structure is important and need to be concruent to obsS:
%           modS.DEM --> the corresponding data within DEM is plotted in a seperate
%           subplot (each DEM goes in a subplot); modS.DEM.data ---> is
%           also a structure with all the datafields in a structure array
%           (at least field defined by fieldname must be included) --> the
%           different data structure arrays are plotted as seperate qq plot
%           in the corresponding subplot (important is that the data
%           structure arrays of the same name are stored in all DEM).
% OPTIONAL INPUT PARAMETERS:
%   
% OUTPUTS:   
%   
% MODIFICATION HISTORY:
%  Thomas Rölli, thomasroelli@gmail.com, 20180612, init

global climada_global
if ~climada_init_vars, return; end

%check arguments
if ~exist('modS') modS = []; end
if ~exist('obsS') obsS = []; end
if ~exist('fieldname') fieldname = []; end
if ~exist('qqOrHist') qqOrHist = []; end
if ~exist('limQQ') limQQ = []; end
if ~exist('binWidth') binWidth = []; end
if ~exist('limHIST') limHIST = []; end
if ~exist('CDF_plot') CDF_plot = []; end
if ~exist('limCDF') limCDF = []; end


if isempty(modS); return; end
if isempty(obsS); return; end
if isempty(fieldname); fieldname = 'length'; end
if isempty(qqOrHist); qqOrHist = 'both'; end 
if isempty(limQQ); limQQ = []; end %is defined below
if isempty(binWidth); binWidth = []; end %is defined below
if isempty(limHIST); limHIST = []; end %is defined below
if isempty(CDF_plot); CDF_plot = 0; end
if isempty(limCDF); limCDF = []; end %is defined below



names = fieldnames(modS);
n_runs = numel(names);


%%%%%qqplots%%%%%
if strcmp(qqOrHist,'both') || strcmp(qqOrHist,'qq')
fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1]);
subplot1(floor((n_runs-1)/3)+1,3,'Gap',[0.02 0.02],'XTickL','Margin','YTickL','Margin','FontS',15)
for i=1:n_runs
    %select corresponding structure
    model = modS.(char(names(i)));
    obser = obsS.(char(names(i)));
    %names of subqq-plots
    names_subqq = fieldnames(model);
    %select corresponding values
    model_val = [model.(char(names_subqq(1))).(fieldname)];
    obser_val = [obser.(char(names_subqq(1))).(fieldname)];
    
    %extract data from qqplot --> to make own plot
    q_fig = figure();
    qq = qqplot(obser_val,model_val);
    qq_obs = qq(1).XData;
    qq_mod = qq(1).YData;
    close(q_fig)
    
    %plot qq
    subplot1(i)
    h(1) = plot(qq_obs,qq_mod,'o','MarkerSize',5,'MarkerEdgeColor','black',...
        'MarkerFaceColor',[0.5 0.5 0.5],'DisplayName',char(names_subqq(1)));
    %if i==1, legend; end
    axis equal
    grid on
    %text as title of plot --> takes fieldnames
    text(0.90,0.95,char(names(i)),'FontSize',15,'FontWeight','bold',...
        'Units','normalized','HorizontalAlignment','right')
    %set limite of axes if not defined yet
    xl = xlim;
    yl = ylim;
    if isempty(limQQ)
       limQQ(1) = min(xl(1),yl(1));
       limQQ(2) = max(xl(2),yl(2));
    end
    xlim(limQQ)
    ylim(limQQ)
    xl = xlim;
    yl = ylim;
    %label of axis just at boarder subplots
    if floor((i-1)/3)+1==floor((n_runs-1)/3)+1, xlabel(['observed quantiles (' fieldname ')'],'FontWeight','bold'); end
    if mod(i-1,3)==0, ylabel(['modelled quantiles'],'FontWeight','bold'); end
    hold on
    %plot frame
    plot([xl(1),xl(2),xl(2),xl(1),xl(1)],[yl(1),yl(1),yl(2),yl(2),yl(1)],...
    '-','LineWidth',2,'color','black')
    %plot deciles
    plot(quantile(obser_val,[0.1:0.1:0.9]),quantile(model_val,[0.1:0.1:0.9]),'o',...
        'MarkerSize',7,'MarkerEdgeColor','black','MarkerFaceColor','red')
    %plot 1:1 line
    plot(limQQ,limQQ,'--','LineWidth',1.5,'Color','black')
    
    markertype = 'x';
    col = [60 120 216;204 65 36;255 153 0;55 118 29;83 193 176;163 58 203]/255;
    for j=2:numel(names_subqq)
        
        %select corresponding values
        model_val = [model.(char(names_subqq(j))).(fieldname)];
        obser_val = [obser.(char(names_subqq(j))).(fieldname)];
        %extract data from qqplot --> to make own plot
        q_fig = figure();
        qq = qqplot(obser_val,model_val);
        qq_obs = qq(1).XData;
        qq_mod = qq(1).YData;
        close(q_fig)

        %plot qq
        %subplot1(i)
        h(j) = plot(qq_obs,qq_mod,markertype,'MarkerSize',5,'MarkerEdgeColor',col(j-1,:),'DisplayName',char(names_subqq(j)));
        uistack(h(j),'bottom');
        axis equal
        grid on
        %text as title of plot --> takes fieldnames
        text(0.90,0.95,char(names(i)),'FontSize',15,'FontWeight','bold',...
            'Units','normalized','HorizontalAlignment','right')
        %set limite of axes if not defined yet
        xl = xlim;
        yl = ylim;
        if isempty(limQQ)
           limQQ(1) = min(xl(1),yl(1));
           limQQ(2) = max(xl(2),yl(2));
        end
        xlim(limQQ)
        ylim(limQQ)
        xl = xlim;
        yl = ylim;
        %label of axis just at boarder subplots
        if floor((i-1)/3)+1==floor((n_runs-1)/3)+1, xlabel(['observed quantiles (' fieldname ')'],'FontWeight','bold'); end
        if mod(i-1,3)==0, ylabel(['modelled quantiles'],'FontWeight','bold'); end
        
    end
    hold off   
end
end

%%%%histograms%%%%
if strcmp(qqOrHist,'both') || strcmp(qqOrHist,'hist')
fig = figure('units','normalized','outerposition',[0 0.3 1 0.65]);
subplot1(floor((n_runs-1)/3)+1,3,'Gap',[0.02 0.02],'XTickL','Margin','YTickL','Margin','FontS',15)
for i=1:n_runs
    %select corresponding structure
    model = modS.(char(names(i)));
    obser = obsS.(char(names(i)));
    %names of subqq-plots
    names_subqq = fieldnames(model);
    %select corresponding values
    model_val = [model.(char(names_subqq(1))).(fieldname)];
    obser_val = [obser.(char(names_subqq(1))).(fieldname)];
    
    %define binWidth if not given
    if isempty(binWidth)
        [~,edges] = histcounts(model_val);
        binWidth = edges(2)-edges(1);
    end
    
    %plots
    subplot1(i)
    histogram(model_val,'BinWidth',binWidth)
    grid on
    %text as title of subplot
    text(0.90,0.95,char(names(i)),'FontSize',15,'FontWeight','bold',...
        'Units','normalized','HorizontalAlignment','right')
    %label of axis just at boarder subplots
    if floor((i-1)/3)+1==floor((n_runs-1)/3)+1, xlabel(fieldname,'FontWeight','bold'); end
    if mod(i-1,3)==0, ylabel('Frequency','FontWeight','bold'); end
    %set limite of axes if not defined yet
    xl = xlim;
    %yl = ylim;
    if isempty(limHIST)
       limHIST = xlim;
    end
    if i==1, yl = ylim; end
    xlim(limHIST)
    ylim(yl)
    xl = xlim;
    yl = ylim;
    %plot observation
    hold on
    histogram(obser_val,'BinWidth',binWidth)
    %plot frame
    plot([xl(1),xl(2),xl(2),xl(1),xl(1)],[yl(1),yl(1),yl(2),yl(2),yl(1)],...
    '-','LineWidth',2,'color','black')

    hold off

end
end

if CDF_plot
fig = figure('units','normalized','outerposition',[0 0.3 1 0.65]);
subplot1(floor((n_runs-1)/3)+1,n_runs,'Gap',[0.02 0.02],'XTickL','Margin','YTickL','Margin','FontS',15)
for i=1:n_runs
    %select corresponding structure
    model = modS.(char(names(i)));
    %names of subqq-plots
    names_dat = fieldnames(model);
    %select corresponding values --> first value
    model_val = [model.(char(names_dat(1))).(fieldname)];

    col = [60 120 216;204 65 36;255 153 0;55 118 29;83 193 176;163 58 203]/255;
    [f_mod,x_mod] = ecdf(model_val);

    %plots
    subplot1(i)
    plot(x_mod,f_mod,'Color',col(1,:),'LineWidth',2,'DisplayName',char(names_dat(1)))
    grid on
    hold on
    xl = xlim;

    %set xlim if not given
    if isempty(limCDF)
       limCDF = xlim;
    end
    if i==1, yl = ylim; end
 
    %iteration through other data --> plot in subplot
    for j=2:numel(names_subqq)
        model_val = [model.(char(names_dat(j))).(fieldname)];
        
        [f_mod,x_mod] = ecdf(model_val);
       
        %plots
        subplot1(i)
        plot(x_mod,f_mod,'Color',col(j,:),'LineWidth',2,'DisplayName',char(names_dat(j)))

    end
    xlim(limCDF)
    ylim(yl)
    xl = xlim;
    yl = ylim;
    
    text(0.90,0.95,char(names(i)),'FontSize',15,'FontWeight','bold',...
        'Units','normalized','HorizontalAlignment','right')
    
    plot([xl(1),xl(2),xl(2),xl(1),xl(1)],[yl(1),yl(1),yl(2),yl(2),yl(1)],...
        '-','LineWidth',2,'color','black')
    
    if floor((i-1)/3)+1==floor((n_runs-1)/3)+1, xlabel([fieldname '[]'],'FontWeight','bold'); end
    if mod(i-1,3)==0, ylabel(['F(' fieldname ')'],'FontWeight','bold'); end
    hold off
end

end