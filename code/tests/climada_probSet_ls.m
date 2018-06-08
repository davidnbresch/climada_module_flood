function [vmax,phi] = climada_probSet_ls(reachAngle,phi_hwidth,edges_vmax,numrun,binORlin,fig)

% Generates probabilistic parameter sets (vmax, phi)
% MODULE:
%  flood
% NAME:
%  climada_probSet_ls
% PURPOSE:
%    
% CALLING SEQUENCE:
%   
% EXAMPLE:
%   
% INPUTS:
%   edges_vmax: vector with 4 elements. First and last element gives range
%               of possible vmax values. From the first to the second edge
%               in the vector, the probability is linear increasing.
%               Between element 2 and 3, the probability is constant and
%               between 3 and 4, linear decreasing. Outside of the range,
%               (element 1 and 4) the probabiliy is zero
%   edges_phi:  Vector which was derived from [~,edges_phi] =
%               histcounts(...). Defines the edges of the histogram of
%               angle of reach (phi) of landslide inventory.
%   count_phi:  Vector which was derived from [count_phi,~] =
%               histcounts(...). Values should be normalized (use
%               'Normalization', 'probability'). Defines bar height of
%               histogram of angle of reach (phi) of landslide inventory.
% OPTIONAL INPUT PARAMETERS:
%   numrun:     Defines the number of random numbers which are generated to
%               produce the probabilistic paramater set. Default = 10^6
%   binORlin:   string with 'bin' or lin' (default = 'lin'). Defines way how probability
%               curve of histogram of phi is constructed.
%               'bin': probability curve is constructed along bars of
%               histogram --> steep steps between bars
%               'lin': probability curve is cunstructed by connecting
%               center points of bar with linear functions.
% OUTPUTS:   
%   para_set:   Matrix with probabilistic paramerters. para_set(:,1) = vmax
%               para_set(:,2) = phi
% MODIFICATION HISTORY:
% Thomas Rölli, thomasroelli@gmail.com, 20180531, init

global climada_global
if ~climada_init_vars, return; end

if ~exist('reachAngle') return; end
if ~exist('phi_hwidth') return; end
if ~exist('edges_vmax') return; end
if ~exist('numrun') numrun = 10^6; end
if ~exist('binORlin') binORlin = 'lin'; end
if ~exist('fig') fig = 0; end

[count_phi,edges_phi] = histcounts(reachAngle,'BinWidth',phi_hwidth,'Normalization', 'probability');
%normalise bars according to area --> area of bars gets equal to 1
count_phi = count_phi/sum(count_phi*phi_hwidth);
count_phi(isnan(count_phi)) = 0;

numbin_vmax = numel(edges_vmax)-1;
numbin_phi = numel(edges_phi)-1;
width = edges_phi(2)-edges_phi(1);

%generate uniform random number within the ranges of the two histograms
%(vmax and phi)
rand_vmax = edges_vmax(1) + rand([1,numrun])*(edges_vmax(end)-edges_vmax(1));
rand_phi = edges_phi(1) + rand([1,numrun])*(edges_phi(end)-edges_phi(1));

%%
%calculate probability for random vmax
prob_vmax = zeros(size(rand_vmax));

%calculate probability
%out of edges = 0;
idx = rand_vmax<edges_vmax(1) | rand_vmax>edges_vmax(end);
prob_vmax(idx) = 0;
%constant probability between 4 and 8 = 1
idx = rand_vmax>=edges_vmax(2) & rand_vmax<=edges_vmax(3);
prob_vmax(idx) = 1;
%linear increasing at left edge (from 1 up to 4)
idx = rand_vmax>=edges_vmax(1) & rand_vmax<edges_vmax(2);
prob_vmax(idx) = (rand_vmax(idx)-edges_vmax(1))/(edges_vmax(2)-edges_vmax(1));
%linear decreasing at right edge (from 8 up to 11)
idx = rand_vmax>edges_vmax(3) & rand_vmax<=edges_vmax(end);
prob_vmax(idx) = 1-(rand_vmax(idx)-edges_vmax(3))/(edges_vmax(end)-edges_vmax(3));

%calculate area under curve to normalize probability
a1 = ((edges_vmax(2)-edges_vmax(1))*1)/2;
a2 = ((edges_vmax(end)-edges_vmax(3))*1)/2;
a3 = (edges_vmax(3)-edges_vmax(2))*1;

prob_vmax = prob_vmax/(a1+a2+a3);
clear a1 a2 a3 idx

%%
%calculate probability for random phi

%stepwise change of probability --> 'bin'
prob_phi_bin = zeros(size(rand_phi));

%set values out of range = 0
idx = rand_phi<edges_phi(1) | rand_phi>=edges_phi(end);
prob_phi_bin(idx) = 0;

for i=1:numel(count_phi)
    idx = rand_phi>=edges_phi(i) & rand_phi<edges_phi(i+1);
    prob_phi_bin(idx) = count_phi(i);
end
clear idx i


%for linear approach --> linear change of probability --> 'lin'
prob_phi_lin = zeros(size(rand_phi));
edges_phi_lin = [edges_phi(1) edges_phi(1:end-1)+width/2 edges_phi(end)];
count_phi_lin = [0 count_phi 0];

idx = rand_phi<edges_phi_lin(1) | rand_phi>=edges_phi_lin(end);
prob_phi_lin(idx) = 0;
a = 0;
for i=1:numel(count_phi_lin)-1
    idx = rand_phi>=edges_phi_lin(i) & rand_phi<edges_phi_lin(i+1);
    prob_phi_lin(idx) = count_phi_lin(i)+((count_phi_lin(i+1)-count_phi_lin(i))...
        /(edges_phi_lin(i+1)-edges_phi_lin(i)))*(rand_phi(idx)-edges_phi_lin(i));
    %a = a+abs((count_phi(i+1)-count_phi(i))*(edges_phi(i+1)-edges_phi(i))/2);
end
%normalize
%prob_phi = prob_phi/a;
clear idx i

%%
%generate uniform random numbers (rand_wig) and look for each set of parameters if
%probability product (prob_vmax*prob_phi) greater than rand_wig --> if yes
%choose as a parameter set if not remove
if binORlin == 'bin'
    prob_phi = prob_phi_bin;
end
if binORlin == 'lin'
    prob_phi = prob_phi_lin;
end

rand_wig = rand([1,numrun]);

idx = find((prob_vmax.*prob_phi) > rand_wig);

vmax = rand_vmax(idx);
phi = rand_phi(idx);

%%
%plots

if fig
    %sort random numbers --> better for plot
    [~,idx_svmax] = sort(rand_vmax);
    [~,idx_sphi] = sort(rand_phi);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%plot density functions of angle of reach and vmax
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot1(2,2,'Gap',[0.02 0.06],'XTickL','All','YTickL','All','FontS',12)

    %plot density angle of reach
    subplot1(1)
    bar(edges_phi(1:end-1)+phi_hwidth/2,count_phi,1,'FaceColor',[.75 .75 .75]);
    grid on
    ylabel('Density','FontWeight','bold');
    xlabel('PHI [\circ]','FontWeight','bold')
    hold on
    a = area(rand_phi(idx_sphi),prob_phi(idx_sphi),'FaceColor',[0, 0.4470, 0.7410]);
    alpha(a,0.25)
    plot(rand_phi(idx_sphi),prob_phi(idx_sphi),'color',[0, 0.4470, 0.7410],'LineWidth',2);
    %ylim([0,0.25])
    yl = ylim;
    xl = xlim;
    %boarder plot
    plot([xl(1) xl(2) xl(2) xl(1) xl(1)],[yl(1) yl(1) yl(2) yl(2) yl(1)],'-','LineWidth',2,'color','black')
    hold off

    %plot density vmax
    subplot1(2)
    a = area(rand_vmax(idx_svmax),prob_vmax(idx_svmax),'FaceColor',[0.8500 0.3250 0.0980]);
    alpha(a,0.25)
    hold on
    plot(rand_vmax(idx_svmax),prob_vmax(idx_svmax),'color',[0.8500 0.3250 0.0980],'LineWidth',2);
    grid on
    xlabel('VMAX [m/s]','FontWeight','bold')
    yl = ylim;
    xl = xlim;
    %boarder plot
    plot([xl(1) xl(2) xl(2) xl(1) xl(1)],[yl(1) yl(1) yl(2) yl(2) yl(1)],'-','LineWidth',2,'color','black')
    ylim(yl)
    hold off
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %figure contour with density 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot1(3)
    cmap = jet(512);
    cmap = cmap(150:450,:);
    
    %get density
    [ dmap ] = dataDensity( vmax, phi, 257, 257, [1 11 5 60]);
    PHI = ([0:1/256:1])*(max(edges_phi)-min(edges_phi))+min(edges_phi);
    VMAX = ([0:1/256:1])*(max(edges_vmax)-min(edges_vmax))+min(edges_vmax);
    [VMAX,PHI] = meshgrid(VMAX,PHI);
    
    p=pcolor(PHI,VMAX,dmap);
    grid on
    xlim([min(phi),max(phi)]);
    ylim([min(vmax),max(vmax)]);
    xlabel('PHI [\circ]','FontWeight','bold')
    ylabel('VMAX [m/s]','FontWeight','bold')
    colormap(cmap)
    set(p, 'EdgeColor', 'none','facealpha',0.75);
    hold on
    plot(phi,vmax,'.','color','black')
    %boarder plot
    xl = xlim;
    yl = ylim;
    plot([xl(1) xl(2) xl(2) xl(1) xl(1)],[yl(1) yl(1) yl(2) yl(2) yl(1)],'-','LineWidth',2,'color','black')
    xlim(xl)
    ylim(yl)
    hold off
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %figure for 3d contour with density on side
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %define colormap for density
    cmap = jet(512);
    cmap = cmap(150:450,:);
    
    figure('units','normalized','outerposition',[0.5 0.1 0.5 0.5])
    plot3(rand_phi(idx_sphi),zeros(size(rand_phi))+max(rand_vmax),prob_phi(idx_sphi),'color',[0, 0.4470, 0.7410],'LineWidth',2)
    set(gca,'fontsize', 12);
    grid on
    xlim([min(edges_phi),max(edges_phi)])
    ylim([min(edges_vmax),max(edges_vmax)])
    hold on
    plot3(zeros(size(rand_vmax))+max(rand_phi),rand_vmax(idx_svmax),prob_vmax(idx_svmax),'color',[0.8500 0.3250 0.0980],'LineWidth',2)
    xl = xlim;
    yl = ylim;
    zl = zlim;
    xlabel('PHI [\circ]','FontWeight','bold','Rotation',11)
    ylabel('VMAX [m/s]','FontWeight','bold','Rotation',-17)
    zlabel('Density','FontWeight','bold')
    %plot boarder
    [ dmap ] = dataDensity( vmax, phi, 257, 257, [1 11 5 60]);
    PHI = ([0:1/256:1])*(max(edges_phi)-min(edges_phi))+min(edges_phi);
    VMAX = ([0:1/256:1])*(max(edges_vmax)-min(edges_vmax))+min(edges_vmax);
    [VMAX,PHI] = meshgrid(VMAX,PHI);
    p=pcolor(PHI,VMAX,dmap);
    colormap(cmap)
    set(p, 'EdgeColor', 'none','facealpha',0.75);
    plot(phi,vmax,'.','MarkerSize',3,'color','black')
    %plot boarder
    plot3([xl(2) xl(2)],[yl(2) yl(2)],[zl(1) zl(2)],'--','color',[0.75 0.75 0.75],'LineWidth',1)
    %plot3(ones(5)*xl(2),[yl(1) yl(2) yl(2) yl(1) yl(1)],[zl(1) zl(1) zl(2) zl(2) zl(1)],'color','black','LineWidth',1)
    %plot3([xl(1) xl(1) xl(2)],[yl(2) yl(1) yl(1)],zeros(3),'color','black','LineWidth',1)
    zlim(zl)
    hold off
    
end




end

