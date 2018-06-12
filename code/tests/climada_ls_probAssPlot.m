function climada_ls_probAssPlot(modS,obsS,fieldname,qqOrHist,limQQ,binWidth,limHIST)

% Plots distributions (qq-plot and histogram) of probabilistically modelled and
% observed (landslide inventory) variable 
% MODULE:
%  
% NAME:
%  climada_probSet_ls
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


if isempty(modS); return; end
if isempty(obsS); return; end
if isempty(fieldname); fieldname = 'length'; end
if isempty(qqOrHist); qqOrHist = 'both'; end 
if isempty(limQQ); limQQ = []; end %is defined below
if isempty(binWidth); binWidth = []; end %is defined below
if isempty(limHIST); limHIST = []; end %is defined below



names = fieldnames(modS);
n_runs = numel(names);


%%%%%qqplots%%%%%
if strcmp(qqOrHist,'both') || strcmp(qqOrHist,'qq')
fig = figure('units','normalized','outerposition',[0 0 1 1]);
subplot1(floor((n_runs-1)/3)+1,3,'Gap',[0.02 0.02],'XTickL','Margin','YTickL','Margin','FontS',13)
for i=1:n_runs
    %select corresponding structure
    model = modS.(char(names(i)));
    obser = obsS.(char(names(i)));
    %select corresponding values
    model_val = [model.(fieldname)];
    obser_val = [obser.(fieldname)];
    
    %extract data from qqplot --> to make own plot
    q_fig = figure();
    qq = qqplot(obser_val,model_val);
    qq_obs = qq(1).XData;
    qq_mod = qq(1).YData;
    close(q_fig)
    
    %plot qq
    subplot1(i)
    plot(qq_obs,qq_mod,'o','MarkerSize',5,'MarkerEdgeColor','black','MarkerFaceColor',[0.5 0.5 0.5])
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
    hold off
end
end

%%%%histograms%%%%
if strcmp(qqOrHist,'both') || strcmp(qqOrHist,'hist')
fig = figure('units','normalized','outerposition',[0 0.3 1 0.6]);
subplot1(floor((n_runs-1)/3)+1,3,'Gap',[0.02 0.02],'XTickL','Margin','YTickL','Margin','FontS',13)
for i=1:n_runs
    %select corresponding structure
    model = modS.(char(names(i)));
    obser = obsS.(char(names(i)));
    %select corresponding values
    model_val = [model.(fieldname)];
    obser_val = [obser.(fieldname)];
    
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
    yl = ylim;
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

end